"""
Tests for the OMICS adapter experiment node redesign.

Tests verify that the adapter correctly creates:
1. Publication nodes (with study metadata merged in)
2. Experiment nodes from experiments block
3. has_experiment edges from Publication -> Experiment
4. tests_coculture_with edges from Experiment -> OrganismTaxon (coculture only)
5. changes_expression_of edges from Experiment -> Gene
6. Correct edge properties (log2FC, p-value, direction, time_point, time_point_order, etc.)
7. Experiment node properties (name, organism_name, treatment_type, etc.)

Note: Organism nodes are created by the CyanorakNcbi adapter (single source of truth).
The OMICS adapter only references organism IDs in tests_coculture_with edges.
"""

import pytest
import tempfile
import os
import yaml
import pandas as pd
from pathlib import Path

from multiomics_kg.adapters.omics_adapter import OMICSAdapter, MultiOMICSAdapter


@pytest.fixture
def temp_data_dir():
    """Create a temporary directory for test data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_de_data():
    """Sample differential expression data."""
    return pd.DataFrame({
        'Synonym': ['PMM0001', 'PMM0002', 'PMM0003', 'PMM0004'],
        'log2_fold_change': [2.5, -1.8, 0.5, -3.2],
        'adjusted_p_value': [0.001, 0.01, 0.05, 0.0001],
    })


@pytest.fixture
def sample_config(temp_data_dir, sample_de_data):
    """Create sample paperconfig.yaml and data file."""
    # Save sample DE data
    data_file = os.path.join(temp_data_dir, 'de_genes.csv')
    sample_de_data.to_csv(data_file, index=False)

    # Create config with experiments block (new format)
    config = {
        'publication': {
            'papername': 'Test Publication 2024',
            'experiments': {
                'coculture_alteromonas_med4_rnaseq': {
                    'name': 'Test DE Analysis',
                    'organism': 'Prochlorococcus MED4',
                    'treatment_condition': 'Coculture with Alteromonas',
                    'control_condition': 'Axenic',
                    'experimental_context': 'in test growth medium under continuous light',
                    'omics_type': 'RNASEQ',
                    'test_type': 'DESeq2',
                    'treatment_type': 'coculture',
                    'medium': 'Pro99',
                    'temperature': '24C',
                    'light_condition': 'continuous_light',
                    'light_intensity': '',
                    'treatment_organism': 'Alteromonas macleodii',
                    'treatment_taxid': 28108,
                    'treatment_assembly_accession': 'GCF_001077695.1',
                    'table_scope': 'all_detected_genes',
                }
            },
            'supplementary_materials': {
                'supp_table_1': {
                    'type': 'csv',
                    'filename': data_file,
                    'statistical_analyses': [
                        {
                            'id': 'test_de_analysis',
                            'experiment': 'coculture_alteromonas_med4_rnaseq',
                            'timepoint': '24h',
                            'timepoint_hours': 24.0,
                            'name_col': 'Synonym',
                            'logfc_col': 'log2_fold_change',
                            'adjusted_p_value_col': 'adjusted_p_value',
                        }
                    ]
                }
            }
        }
    }

    config_file = os.path.join(temp_data_dir, 'paperconfig.yaml')
    with open(config_file, 'w') as f:
        yaml.dump(config, f)

    return config_file


@pytest.fixture
def adapter_with_mock_extracted_data(sample_config):
    """Create adapter with mocked extracted_data to skip PDF extraction."""
    adapter = OMICSAdapter(config_file=sample_config)
    # Mock the extracted_data that would normally come from PDF extraction
    adapter.extracted_data = {
        'publication': {
            'publication_id': 'test_pub_2024',
            'title': 'Test Publication',
            'doi': '10.1234/test.2024',
            'authors': ['Test Author'],
            'description': 'Test study description',
            'study_type': 'Transcriptomics',
            'organism': ['Prochlorococcus MED4'],
        },
    }
    return adapter


class TestPublicationNodeCreation:
    """Test that publication nodes are created correctly."""

    def test_publication_node_created(self, adapter_with_mock_extracted_data):
        """Verify publication node is created with merged study fields."""
        adapter = adapter_with_mock_extracted_data
        nodes = adapter.get_nodes()

        pub_nodes = [n for n in nodes if n[1] == 'publication']
        assert len(pub_nodes) == 1, "Expected exactly one publication node"

        node_id, label, properties = pub_nodes[0]
        assert properties.get('title') == 'Test Publication'
        assert properties.get('description') == 'Test study description'
        assert properties.get('study_type') == 'Transcriptomics'
        # organism is no longer set by the adapter; it's computed post-import as 'organisms'
        assert 'organism' not in properties


class TestOrganismNodeCreation:
    """Test that OMICS adapter does NOT create organism nodes.

    Organism nodes are created by the CyanorakNcbi adapter (single source of truth).
    """

    def test_no_organism_nodes_from_omics(self, adapter_with_mock_extracted_data):
        """Verify OMICS adapter does not create organism nodes."""
        adapter = adapter_with_mock_extracted_data
        nodes = adapter.get_nodes()

        organism_nodes = [n for n in nodes if n[1] == 'organism']
        assert len(organism_nodes) == 0, \
            "OMICS adapter should not create organism nodes (CyanorakNcbi is the single source)"

    def test_no_organism_nodes_without_treatment_taxid(self, temp_data_dir, sample_de_data):
        """Verify no organism node is created regardless of config."""
        data_file = os.path.join(temp_data_dir, 'de_genes.csv')
        sample_de_data.to_csv(data_file, index=False)

        config = {
            'publication': {
                'papername': 'Test Publication 2024',
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [
                            {
                                'type': 'RNASEQ',
                                'name': 'Test DE Analysis',
                                'id': 'test_de_analysis',
                                'test_type': 'DESeq2',
                                'organism': 'Prochlorococcus MED4',
                                'name_col': 'Synonym',
                                'logfc_col': 'log2_fold_change',
                                'adjusted_p_value_col': 'adjusted_p_value',
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_no_treatment.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {'publication_id': 'test', 'doi': '10.1234/test'},
        }

        nodes = adapter.get_nodes()
        organism_nodes = [n for n in nodes if n[1] == 'organism']
        assert len(organism_nodes) == 0


class TestExperimentNodeCreation:
    """Test that experiment nodes are created correctly."""

    def test_experiment_node_created(self, adapter_with_mock_extracted_data):
        """Verify experiment node is created from experiments block."""
        adapter = adapter_with_mock_extracted_data
        nodes = adapter.get_nodes()

        exp_nodes = [n for n in nodes if n[1] == 'experiment']

        assert len(exp_nodes) == 1, f"Expected exactly one experiment node, got {len(exp_nodes)}"

        node_id, label, properties = exp_nodes[0]

        # Check that ID includes publication ID and experiment key
        assert 'test_pub_2024' in node_id or '10.1234' in node_id, \
            f"Expected experiment node ID to include publication ID, got {node_id}"
        assert 'coculture_alteromonas_med4_rnaseq' in node_id, \
            f"Expected experiment node ID to include experiment key, got {node_id}"

        # Check properties
        assert properties.get('name') == 'Test DE Analysis', \
            f"Expected name 'Test DE Analysis', got {properties.get('name')}"
        assert properties.get('organism_name') == 'Prochlorococcus MED4', \
            f"Expected organism_name 'Prochlorococcus MED4', got {properties.get('organism_name')}"
        assert properties.get('treatment_type') == ['coculture'], \
            f"Expected treatment_type ['coculture'], got {properties.get('treatment_type')}"
        assert properties.get('light_condition') == 'continuous_light', \
            f"Expected light_condition 'continuous_light', got {properties.get('light_condition')}"
        assert properties.get('temperature') == '24C', \
            f"Expected temperature '24C', got {properties.get('temperature')}"

    def test_no_experiment_node_without_experiments_block(self, temp_data_dir, sample_de_data):
        """Verify no experiment node is created when experiments block is missing."""
        data_file = os.path.join(temp_data_dir, 'de_genes.csv')
        sample_de_data.to_csv(data_file, index=False)

        config = {
            'publication': {
                'papername': 'Test Publication 2024',
                # No experiments block
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [
                            {
                                'id': 'test_de_analysis',
                                'name_col': 'Synonym',
                                'logfc_col': 'log2_fold_change',
                                'adjusted_p_value_col': 'adjusted_p_value',
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_no_exp.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {'publication_id': 'test', 'doi': '10.1234/test'},
        }

        nodes = adapter.get_nodes()
        exp_nodes = [n for n in nodes if n[1] == 'experiment']

        assert len(exp_nodes) == 0, "No experiment node should be created without experiments block"


class TestExperimentToGeneEdges:
    """Test that experiment → gene expression edges are created correctly.

    The sample_config fixture uses the new experiments block format. All expression
    edges should have label 'changes_expression_of' with experiment as source.
    """

    def test_expression_edge_created(self, adapter_with_mock_extracted_data):
        """Verify changes_expression_of edges are created."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        assert len(expression_edges) == 4, f"Expected 4 expression edges (one per gene), got {len(expression_edges)}"

    def test_edge_source_is_experiment(self, adapter_with_mock_extracted_data):
        """Verify edge source is the experiment node."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            # Source should be the experiment ID (pubid_experiment_key)
            assert 'coculture_alteromonas_med4_rnaseq' in source_id, \
                f"Expected source to contain experiment key, got {source_id}"

    def test_edge_target_is_gene(self, adapter_with_mock_extracted_data):
        """Verify edge target is a gene."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'ncbigene' in target_id.lower() or 'PMM' in target_id, \
                f"Expected target to be gene ID, got {target_id}"

    def test_edge_has_log2_fold_change(self, adapter_with_mock_extracted_data):
        """Verify edges have log2_fold_change property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'log2_fold_change' in properties, \
                f"Expected log2_fold_change in edge properties, got {properties.keys()}"
            assert isinstance(properties['log2_fold_change'], float), \
                f"Expected log2_fold_change to be float, got {type(properties['log2_fold_change'])}"

    def test_edge_has_adjusted_p_value(self, adapter_with_mock_extracted_data):
        """Verify edges have adjusted_p_value property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'adjusted_p_value' in properties, \
                f"Expected adjusted_p_value in edge properties, got {properties.keys()}"

    def test_edge_has_expression_direction(self, adapter_with_mock_extracted_data):
        """Verify edges have expression_direction property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'expression_direction' in properties, \
                f"Expected expression_direction in edge properties"
            assert properties['expression_direction'] in ['up', 'down'], \
                f"Expected expression_direction to be 'up' or 'down', got {properties['expression_direction']}"

    def test_expression_direction_matches_fold_change(self, adapter_with_mock_extracted_data):
        """Verify expression_direction matches sign of log2_fold_change."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            fc = properties['log2_fold_change']
            direction = properties['expression_direction']

            if fc > 0:
                assert direction == 'up', f"Positive fold change ({fc}) should have direction 'up', got '{direction}'"
            else:
                assert direction == 'down', f"Negative fold change ({fc}) should have direction 'down', got '{direction}'"

    def test_edge_has_time_point(self, adapter_with_mock_extracted_data):
        """Verify edges have time_point property when specified in config."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'time_point' in properties, \
                f"Expected time_point in edge properties"
            assert properties['time_point'] == '24h', \
                f"Expected time_point '24h', got {properties['time_point']}"

    def test_edge_has_time_point_order(self, adapter_with_mock_extracted_data):
        """Verify edges have time_point_order property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'time_point_order' in properties, \
                f"Expected time_point_order in edge properties"
            assert properties['time_point_order'] == 1, \
                f"Expected time_point_order 1, got {properties['time_point_order']}"

    def test_edge_has_time_point_hours(self, adapter_with_mock_extracted_data):
        """Verify edges have time_point_hours property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'time_point_hours' in properties, \
                f"Expected time_point_hours in edge properties"
            assert properties['time_point_hours'] == 24.0, \
                f"Expected time_point_hours 24.0, got {properties['time_point_hours']}"

    def test_edge_no_old_metadata(self, adapter_with_mock_extracted_data):
        """Verify expression edges do not carry old metadata (now on experiment node)."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            # These should NOT be on the edge anymore (they're on the experiment node)
            assert 'publications' not in properties
            assert 'omics_type' not in properties
            assert 'organism_name' not in properties
            assert 'treatment_condition' not in properties
            assert 'control_condition' not in properties
            assert 'experimental_context' not in properties
            assert 'statistical_test' not in properties
            assert 'analysis_name' not in properties


class TestNonCocultureExperimentEdges:
    """Test expression edges for non-coculture experiments (no treatment_organism)."""

    @pytest.fixture
    def env_experiment_config(self, temp_data_dir, sample_de_data):
        """Create config with a non-coculture experiment (gas shock)."""
        data_file = os.path.join(temp_data_dir, 'de_genes.csv')
        sample_de_data.to_csv(data_file, index=False)

        config = {
            'publication': {
                'papername': 'Test Env Experiment Publication 2024',
                'experiments': {
                    'gas_shock_air': {
                        'name': 'Gene expression response to air',
                        'organism': 'Prochlorococcus MED4',
                        'treatment_condition': '0.036% CO2 / 21% O2 (air)',
                        'control_condition': 'Time 0 (pre-shock)',
                        'experimental_context': 'Pro99 medium, constant light',
                        'omics_type': 'MICROARRAY',
                        'test_type': 'Affymetrix microarray',
                        'treatment_type': 'gas_shock',
                        'medium': 'Pro99',
                        'temperature': '24C',
                        'light_condition': 'continuous light',
                        'light_intensity': '',
                    }
                },
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [
                            {
                                'id': 'gas_shock_air_analysis',
                                'experiment': 'gas_shock_air',
                                'name_col': 'Synonym',
                                'logfc_col': 'log2_fold_change',
                                'adjusted_p_value_col': 'adjusted_p_value',
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_env_exp.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        return config_file

    @pytest.fixture
    def env_adapter(self, env_experiment_config):
        """Create adapter with non-coculture experiment config."""
        adapter = OMICSAdapter(config_file=env_experiment_config)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'test_env_pub_2024',
                'title': 'Test Env Publication',
                'doi': '10.1234/test_env.2024',
            },
        }
        return adapter

    def test_expression_edges_created(self, env_adapter):
        """Verify changes_expression_of edges are created for non-coculture experiment."""
        edges = env_adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        assert len(expression_edges) == 4, f"Expected 4 expression edges, got {len(expression_edges)}"

    def test_edge_source_is_experiment(self, env_adapter):
        """Verify edge source is the experiment node ID."""
        edges = env_adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'gas_shock_air' in source_id, \
                f"Expected source to contain experiment key, got {source_id}"

    def test_edge_properties(self, env_adapter):
        """Verify edges have correct properties."""
        edges = env_adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'log2_fold_change' in properties
            assert 'adjusted_p_value' in properties
            assert 'expression_direction' in properties
            assert properties['expression_direction'] in ['up', 'down']

    def test_no_tests_coculture_with_edge(self, env_adapter):
        """Verify no tests_coculture_with edge when experiment has no treatment_organism."""
        edges = env_adapter.get_edges()
        coculture_edges = [e for e in edges if e[3] == 'tests_coculture_with']
        assert len(coculture_edges) == 0, \
            f"Expected no tests_coculture_with edges for non-coculture experiment, got {len(coculture_edges)}"


class TestSkipRowsWithMissingFoldChange:
    """Test that rows with missing/empty fold change values are skipped."""

    @pytest.fixture
    def de_data_with_missing_fc(self):
        """Sample data with missing, empty, and NA fold change values."""
        return pd.DataFrame({
            'Synonym': ['PMM0001', 'PMM0002', 'PMM0003', 'PMM0004', 'PMM0005'],
            'log2_fold_change': [2.5, '', float('nan'), 'NA', -1.0],
            'adjusted_p_value': [0.001, 0.01, 0.05, 0.0001, 0.02],
        })

    @pytest.fixture
    def adapter_with_missing_fc(self, temp_data_dir, de_data_with_missing_fc):
        """Create adapter with data containing missing fold change values."""
        data_file = os.path.join(temp_data_dir, 'de_genes_missing_fc.csv')
        de_data_with_missing_fc.to_csv(data_file, index=False)

        config = {
            'publication': {
                'papername': 'Test Missing FC 2024',
                'experiments': {
                    'test_exp': {
                        'name': 'Test DE Analysis',
                        'organism': 'Prochlorococcus MED4',
                        'treatment_condition': 'Treatment',
                        'control_condition': 'Control',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'coculture',
                        'treatment_organism': 'Alteromonas macleodii',
                        'treatment_taxid': 28108,
                        'treatment_assembly_accession': 'GCF_001077695.1',
                    }
                },
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [
                            {
                                'id': 'test_missing_fc',
                                'experiment': 'test_exp',
                                'name_col': 'Synonym',
                                'logfc_col': 'log2_fold_change',
                                'adjusted_p_value_col': 'adjusted_p_value',
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_missing_fc.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'test_missing_fc',
                'title': 'Test',
                'doi': '10.1234/test_fc',
            },
        }
        return adapter

    def test_only_valid_fc_rows_create_edges(self, adapter_with_missing_fc):
        """Verify that only rows with valid fold change values create edges."""
        edges = adapter_with_missing_fc.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        # Only PMM0001 (2.5) and PMM0005 (-1.0) have valid FC values
        assert len(expression_edges) == 2, \
            f"Expected 2 edges (skipping empty, NaN, NA), got {len(expression_edges)}"

    def test_skipped_rows_do_not_appear_in_edges(self, adapter_with_missing_fc):
        """Verify that genes with missing FC are not in the edge targets."""
        edges = adapter_with_missing_fc.get_edges()
        target_ids = [e[2] for e in edges if e[3] == 'changes_expression_of']
        target_str = ' '.join(target_ids)
        assert 'PMM0002' not in target_str, "PMM0002 (empty FC) should be skipped"
        assert 'PMM0003' not in target_str, "PMM0003 (NaN FC) should be skipped"
        assert 'PMM0004' not in target_str, "PMM0004 (NA FC) should be skipped"


class TestSkipRows:
    """Test that skip_rows option correctly skips header rows when loading CSV."""

    @pytest.fixture
    def csv_with_multi_row_header(self, temp_data_dir):
        """Create a CSV file with multiple header rows (like Biller 2016 format)."""
        csv_content = (
            "Table S2. Description line\n"
            ",,\n"
            ",,Time past addition:\n"
            "Gene_ID,log2FC,Description\n"
            "GENE_A,-0.212,chaperone\n"
            "GENE_B,0.425,oxidoreductase\n"
            "GENE_C,-0.37,hypothetical\n"
        )
        data_file = os.path.join(temp_data_dir, 'multi_header.csv')
        with open(data_file, 'w') as f:
            f.write(csv_content)
        return data_file

    @pytest.fixture
    def adapter_skip_rows(self, temp_data_dir, csv_with_multi_row_header):
        """Create adapter with skip_rows=3 to skip multi-row header."""
        config = {
            'publication': {
                'papername': 'Test Skip Rows 2024',
                'experiments': {
                    'skip_rows_exp': {
                        'name': 'Test DE with skip_rows',
                        'organism': 'Prochlorococcus NATL2A',
                        'treatment_condition': 'Coculture',
                        'control_condition': 'Axenic',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'coculture',
                        'treatment_organism': 'Alteromonas macleodii',
                        'treatment_taxid': 28108,
                        'treatment_assembly_accession': 'GCF_001077695.1',
                    }
                },
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': csv_with_multi_row_header,
                        'statistical_analyses': [
                            {
                                'id': 'test_skip_rows',
                                'experiment': 'skip_rows_exp',
                                'name_col': 'Gene_ID',
                                'logfc_col': 'log2FC',
                                'adjusted_p_value_col': None,
                                'skip_rows': 3,
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_skip_rows.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'test_skip_rows',
                'title': 'Test Skip Rows',
                'doi': '10.1234/test_skip',
            },
        }
        return adapter

    def test_skip_rows_loads_correct_data(self, adapter_skip_rows):
        """Verify that skip_rows correctly skips multi-row headers."""
        edges = adapter_skip_rows.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        assert len(expression_edges) == 3, \
            f"Expected 3 edges (GENE_A, GENE_B, GENE_C), got {len(expression_edges)}"

    def test_skip_rows_correct_gene_ids(self, adapter_skip_rows):
        """Verify the correct gene IDs are parsed after skipping rows."""
        edges = adapter_skip_rows.get_edges()
        target_ids = [e[2] for e in edges if e[3] == 'changes_expression_of']
        target_str = ' '.join(target_ids)
        assert 'GENE_A' in target_str
        assert 'GENE_B' in target_str
        assert 'GENE_C' in target_str

    def test_skip_rows_correct_fold_change(self, adapter_skip_rows):
        """Verify fold change values are correctly parsed after skipping rows."""
        edges = adapter_skip_rows.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        fc_values = sorted([e[4]['log2_fold_change'] for e in expression_edges])
        assert fc_values == pytest.approx(sorted([-0.37, -0.212, 0.425]))

    def test_no_skip_rows_fails_gracefully(self, temp_data_dir, csv_with_multi_row_header):
        """Verify that without skip_rows, multi-row header CSV fails to find columns."""
        config = {
            'publication': {
                'papername': 'Test No Skip 2024',
                'experiments': {
                    'no_skip_exp': {
                        'name': 'Test No Skip',
                        'organism': 'Prochlorococcus MED4',
                        'treatment_condition': 'Treatment',
                        'control_condition': 'Control',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'coculture',
                        'treatment_organism': 'Alteromonas macleodii',
                        'treatment_taxid': 28108,
                        'treatment_assembly_accession': 'GCF_001077695.1',
                    }
                },
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': csv_with_multi_row_header,
                        'statistical_analyses': [
                            {
                                'id': 'test_no_skip_rows',
                                'experiment': 'no_skip_exp',
                                'name_col': 'Gene_ID',
                                'logfc_col': 'log2FC',
                                'adjusted_p_value_col': None,
                                # No skip_rows — should not find Gene_ID column
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_no_skip.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'test_no_skip',
                'title': 'Test',
                'doi': '10.1234/no_skip',
            },
        }
        edges = adapter.get_edges()
        # Without skip_rows, the column Gene_ID won't be found → no expression edges
        # (A published_expression_data_about edge may still be emitted for the source node)
        expr_edges = [
            e for e in edges
            if e[3] == 'changes_expression_of'
        ]
        assert len(expr_edges) == 0


class TestPvalueAsteriskInLogfc:
    """Test that pvalue_asterisk_in_logfc correctly parses asterisk-marked significance."""

    @pytest.fixture
    def csv_with_asterisks(self, temp_data_dir):
        """Create a CSV with asterisks in logFC column indicating significance."""
        data = pd.DataFrame({
            'Gene_ID': ['GENE_A', 'GENE_B', 'GENE_C', 'GENE_D'],
            'log2FC': ['-0.212*', '0.425', '-0.37*', '0.15'],
        })
        data_file = os.path.join(temp_data_dir, 'asterisk_data.csv')
        data.to_csv(data_file, index=False)
        return data_file

    @pytest.fixture
    def adapter_asterisk(self, temp_data_dir, csv_with_asterisks):
        """Create adapter with pvalue_asterisk_in_logfc=true."""
        config = {
            'publication': {
                'papername': 'Test Asterisk 2024',
                'experiments': {
                    'asterisk_exp': {
                        'name': 'Test DE with asterisk pvalue',
                        'organism': 'Prochlorococcus NATL2A',
                        'treatment_condition': 'Coculture',
                        'control_condition': 'Axenic',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'coculture',
                        'treatment_organism': 'Alteromonas macleodii',
                        'treatment_taxid': 28108,
                        'treatment_assembly_accession': 'GCF_001077695.1',
                    }
                },
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': csv_with_asterisks,
                        'statistical_analyses': [
                            {
                                'id': 'test_asterisk',
                                'experiment': 'asterisk_exp',
                                'name_col': 'Gene_ID',
                                'logfc_col': 'log2FC',
                                'adjusted_p_value_col': None,
                                'pvalue_asterisk_in_logfc': True,
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_asterisk.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'test_asterisk',
                'title': 'Test Asterisk',
                'doi': '10.1234/test_asterisk',
            },
        }
        return adapter

    def test_all_rows_create_edges(self, adapter_asterisk):
        """Verify all rows create edges (asterisk is stripped, not treated as invalid)."""
        edges = adapter_asterisk.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        assert len(expression_edges) == 4, \
            f"Expected 4 edges, got {len(expression_edges)}"

    def test_fold_change_values_parsed_correctly(self, adapter_asterisk):
        """Verify asterisks are stripped and fold change is parsed as float."""
        edges = adapter_asterisk.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        fc_values = sorted([e[4]['log2_fold_change'] for e in expression_edges])
        assert fc_values == pytest.approx(sorted([-0.212, 0.425, -0.37, 0.15]))

    def test_significant_rows_get_low_pvalue(self, adapter_asterisk):
        """Verify asterisk-marked rows get adjusted_p_value < 0.5."""
        edges = adapter_asterisk.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, _, target_id, _, properties = edge
            if 'GENE_A' in target_id or 'GENE_C' in target_id:
                assert 'adjusted_p_value' in properties
                assert properties['adjusted_p_value'] < 0.5, \
                    f"Significant gene should have low p-value, got {properties['adjusted_p_value']}"

    def test_non_significant_rows_get_high_pvalue(self, adapter_asterisk):
        """Verify non-asterisk rows get adjusted_p_value = 1.0."""
        edges = adapter_asterisk.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, _, target_id, _, properties = edge
            if 'GENE_B' in target_id or 'GENE_D' in target_id:
                assert 'adjusted_p_value' in properties
                assert properties['adjusted_p_value'] == 1.0, \
                    f"Non-significant gene should have p-value 1.0, got {properties['adjusted_p_value']}"

    def test_expression_direction_correct(self, adapter_asterisk):
        """Verify expression direction is correct after asterisk stripping."""
        edges = adapter_asterisk.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, _, target_id, _, properties = edge
            fc = properties['log2_fold_change']
            direction = properties['expression_direction']
            if fc > 0:
                assert direction == 'up'
            else:
                assert direction == 'down'

    def test_without_asterisk_flag_stars_still_stripped(self, temp_data_dir, csv_with_asterisks):
        """Verify that without pvalue_asterisk_in_logfc, stars are still stripped from logFC
        so all rows parse successfully (no significance marking via asterisk though)."""
        config = {
            'publication': {
                'papername': 'Test No Asterisk Flag 2024',
                'experiments': {
                    'no_ast_exp': {
                        'name': 'Test No Asterisk',
                        'organism': 'Prochlorococcus MED4',
                        'treatment_condition': 'Treatment',
                        'control_condition': 'Control',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'coculture',
                        'treatment_organism': 'Alteromonas macleodii',
                        'treatment_taxid': 28108,
                        'treatment_assembly_accession': 'GCF_001077695.1',
                    }
                },
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': csv_with_asterisks,
                        'statistical_analyses': [
                            {
                                'id': 'test_no_asterisk',
                                'experiment': 'no_ast_exp',
                                'name_col': 'Gene_ID',
                                'logfc_col': 'log2FC',
                                'adjusted_p_value_col': None,
                                # No pvalue_asterisk_in_logfc — stars are stripped but not
                                # used for significance
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_no_asterisk.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'test_no_ast',
                'title': 'Test',
                'doi': '10.1234/no_ast',
            },
        }
        edges = adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        # All 4 rows should now produce edges — stars are always stripped before float parse
        assert len(expression_edges) == 4, \
            f"Expected 4 edges (stars stripped even without pvalue_asterisk_in_logfc flag), got {len(expression_edges)}"
        # Without the flag, no asterisk-based p-value should be set
        for edge in expression_edges:
            props = edge[4]
            assert 'adjusted_p_value' not in props, \
                "Should not have adjusted_p_value when pvalue_asterisk_in_logfc is False"

    def test_leading_asterisk_parsed_correctly(self, temp_data_dir):
        """Verify leading-asterisk format '* 1.1' is handled correctly."""
        data = pd.DataFrame({
            'Gene_ID': ['GENE_A', 'GENE_B', 'GENE_C', 'GENE_D'],
            'log2FC': ['* -0.212', '0.425', '* -0.37', '0.15'],
        })
        data_file = os.path.join(temp_data_dir, 'leading_asterisk.csv')
        data.to_csv(data_file, index=False)

        config = {
            'publication': {
                'papername': 'Test Leading Asterisk 2024',
                'experiments': {
                    'lead_ast_exp': {
                        'name': 'Test Leading Asterisk',
                        'organism': 'Prochlorococcus MED4',
                        'treatment_condition': 'Treatment',
                        'control_condition': 'Control',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'coculture',
                        'treatment_organism': 'Alteromonas macleodii',
                        'treatment_taxid': 28108,
                        'treatment_assembly_accession': 'GCF_001077695.1',
                    }
                },
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [
                            {
                                'id': 'test_leading_asterisk',
                                'experiment': 'lead_ast_exp',
                                'name_col': 'Gene_ID',
                                'logfc_col': 'log2FC',
                                'adjusted_p_value_col': None,
                                'pvalue_asterisk_in_logfc': True,
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_leading_ast.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'test_lead_ast',
                'title': 'Test',
                'doi': '10.1234/lead_ast',
            },
        }
        edges = adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        # All 4 rows should produce edges
        assert len(expression_edges) == 4, \
            f"Expected 4 edges from leading-asterisk data, got {len(expression_edges)}"

        # Fold change values should be parsed correctly (no asterisk residue)
        fc_values = sorted([e[4]['log2_fold_change'] for e in expression_edges])
        assert fc_values == pytest.approx(sorted([-0.212, 0.425, -0.37, 0.15]))

        # GENE_A and GENE_C (leading *) should be significant; GENE_B and GENE_D not
        for edge in expression_edges:
            _, _, target_id, _, props = edge
            if 'GENE_A' in target_id or 'GENE_C' in target_id:
                assert props.get('adjusted_p_value', 1.0) < 0.5, \
                    f"Leading-asterisk gene should be significant, got {props.get('adjusted_p_value')}"
            if 'GENE_B' in target_id or 'GENE_D' in target_id:
                assert props.get('adjusted_p_value') == 1.0, \
                    f"Non-asterisk gene should have p-value 1.0, got {props.get('adjusted_p_value')}"


class TestSkipRowsWithAsterisk:
    """Test skip_rows and pvalue_asterisk_in_logfc combined (like Biller 2016 supp_table_2)."""

    @pytest.fixture
    def biller_style_csv(self, temp_data_dir):
        """Create a CSV mimicking Biller 2016 format: multi-row header + asterisk significance."""
        csv_content = (
            "Table S2. All genes with significantly different transcript abundances\n"
            ",,\n"
            ",,Time past addition:\n"
            "Gene_ID,24 hours,Description\n"
            "PMN2A_1344,-0.212*,chaperone DnaJ\n"
            "PMN2A_1346,-0.425*,ribosome GTPase\n"
            "PMN2A_1362,0.048,cobalt methyltransferase\n"
            "PMN2A_1375,-0.37*,flavin oxidoreductase\n"
        )
        data_file = os.path.join(temp_data_dir, 'biller_style.csv')
        with open(data_file, 'w') as f:
            f.write(csv_content)
        return data_file

    @pytest.fixture
    def adapter_combined(self, temp_data_dir, biller_style_csv):
        """Create adapter with both skip_rows and pvalue_asterisk_in_logfc."""
        config = {
            'publication': {
                'papername': 'Biller 2016',
                'experiments': {
                    'biller_coculture_exp': {
                        'name': 'DE coculture vs axenic',
                        'organism': 'Prochlorococcus NATL2A',
                        'treatment_condition': 'Coculture',
                        'control_condition': 'Axenic',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'coculture',
                        'treatment_organism': 'Alteromonas macleodii',
                        'treatment_taxid': 28108,
                        'treatment_assembly_accession': 'GCF_001077695.1',
                    }
                },
                'supplementary_materials': {
                    'supp_table_2': {
                        'type': 'csv',
                        'filename': biller_style_csv,
                        'statistical_analyses': [
                            {
                                'id': 'de_biller_style',
                                'experiment': 'biller_coculture_exp',
                                'name_col': 'Gene_ID',
                                'logfc_col': '24 hours',
                                'adjusted_p_value_col': None,
                                'skip_rows': 3,
                                'pvalue_asterisk_in_logfc': True,
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_combined.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'biller_2016_test',
                'title': 'Biller 2016 Test',
                'doi': '10.1234/biller',
            },
        }
        return adapter

    def test_combined_creates_all_edges(self, adapter_combined):
        """Verify all 4 rows produce edges with combined skip_rows + asterisk."""
        edges = adapter_combined.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        assert len(expression_edges) == 4

    def test_combined_correct_fold_changes(self, adapter_combined):
        """Verify fold change values are correct with combined features."""
        edges = adapter_combined.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']
        fc_values = sorted([e[4]['log2_fold_change'] for e in expression_edges])
        assert fc_values == pytest.approx(sorted([-0.425, -0.37, -0.212, 0.048]))

    def test_combined_significance_markers(self, adapter_combined):
        """Verify significant vs non-significant p-values with combined features."""
        edges = adapter_combined.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        significant_count = sum(
            1 for e in expression_edges if e[4].get('adjusted_p_value', 1.0) < 0.5
        )
        non_significant_count = sum(
            1 for e in expression_edges if e[4].get('adjusted_p_value', 1.0) >= 0.5
        )
        assert significant_count == 3, f"Expected 3 significant edges, got {significant_count}"
        assert non_significant_count == 1, f"Expected 1 non-significant edge, got {non_significant_count}"

    def test_combined_edge_source_is_experiment(self, adapter_combined):
        """Verify edge source is experiment node with combined features."""
        edges = adapter_combined.get_edges()
        expression_edges = [e for e in edges if e[3] == 'changes_expression_of']

        for edge in expression_edges:
            _, source_id, _, _, _ = edge
            assert 'biller_coculture_exp' in source_id


class TestEdgeWithRealData:
    """Test with the real Aharonovich 2016 config if available."""

    @pytest.fixture
    def real_config_path(self):
        """Path to real config file."""
        return 'data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml'

    def test_real_config_no_organism_nodes(self, real_config_path):
        """OMICS adapter should not create organism nodes (CyanorakNcbi handles that)."""
        if not os.path.exists(real_config_path):
            pytest.skip("Real config file not found")

        adapter = OMICSAdapter(config_file=real_config_path)
        adapter.download_data(cache=True)

        nodes = adapter.get_nodes()
        organism_nodes = [n for n in nodes if n[1] == 'organism']
        assert len(organism_nodes) == 0, \
            "OMICS adapter should not create organism nodes"

    def test_real_config_creates_experiment_nodes(self, real_config_path):
        """Test that real config creates experiment nodes."""
        if not os.path.exists(real_config_path):
            pytest.skip("Real config file not found")

        adapter = OMICSAdapter(config_file=real_config_path)
        adapter.download_data(cache=True)

        nodes = adapter.get_nodes()
        exp_nodes = [n for n in nodes if n[1] == 'experiment']

        # Aharonovich 2016 has 5 experiments
        assert len(exp_nodes) >= 1, "Expected at least one experiment node"

        # Check that node has expected properties
        exp_node = exp_nodes[0]
        node_id, label, properties = exp_node
        assert 'name' in properties, \
            f"Expected name in experiment node properties, got {properties.keys()}"

    def test_real_config_creates_expression_edges(self, real_config_path):
        """Test that real config creates expression edges."""
        if not os.path.exists(real_config_path):
            pytest.skip("Real config file not found")

        adapter = OMICSAdapter(config_file=real_config_path, test_mode=True)
        adapter.download_data(cache=True)

        edges = adapter.get_edges()
        expression_edges = [
            e for e in edges
            if e[3] == 'changes_expression_of'
        ]

        # Should have expression edges
        assert len(expression_edges) > 0, "Expected expression edges from real config"

        # Verify structure of first edge
        first_edge = expression_edges[0]
        _, source_id, target_id, label, properties = first_edge

        assert label == 'changes_expression_of'
        # Source is experiment ID, not organism
        assert 'coculture_alteromonas_hot1a3_med4_rnaseq' in source_id, \
            f"Expected source to be experiment ID, got {source_id}"
        assert 'log2_fold_change' in properties
        assert 'expression_direction' in properties


class TestMultiOMICSAdapter:
    """Test that MultiOMICSAdapter correctly aggregates results from multiple configs."""

    @pytest.fixture
    def two_config_list(self, temp_data_dir):
        """Create two separate paperconfig.yaml files and a list file pointing to them."""
        configs = []
        for i, (paper, organism, taxid, accession) in enumerate([
            ('Paper A 2024', 'Alteromonas macleodii', 28108, 'GCF_001077695.1'),
            ('Paper B 2024', 'Synechococcus sp.', 32049, 'GCF_000032049.1'),
        ]):
            # Create data file for this config
            data_file = os.path.join(temp_data_dir, f'de_genes_{i}.csv')
            pd.DataFrame({
                'Synonym': [f'GENE_{i}_1', f'GENE_{i}_2'],
                'log2_fold_change': [1.5, -2.0],
                'adjusted_p_value': [0.01, 0.05],
            }).to_csv(data_file, index=False)

            exp_key = f'exp_{organism.replace(" ", "_").replace(".", "")}'
            config = {
                'publication': {
                    'papername': paper,
                    'experiments': {
                        exp_key: {
                            'name': f'Test {paper}',
                            'organism': 'Prochlorococcus MED4',
                            'treatment_condition': 'Treatment',
                            'control_condition': 'Control',
                            'omics_type': 'RNASEQ',
                            'test_type': 'DESeq2',
                            'treatment_type': 'coculture',
                            'treatment_organism': organism,
                            'treatment_taxid': taxid,
                            'treatment_assembly_accession': accession,
                        }
                    },
                    'supplementary_materials': {
                        'supp_table_1': {
                            'type': 'csv',
                            'filename': data_file,
                            'statistical_analyses': [
                                {
                                    'id': f'test_multi_adapter_{organism.replace(" ", "_")}',
                                    'experiment': exp_key,
                                    'name_col': 'Synonym',
                                    'logfc_col': 'log2_fold_change',
                                    'adjusted_p_value_col': 'adjusted_p_value',
                                }
                            ]
                        }
                    }
                }
            }

            config_file = os.path.join(temp_data_dir, f'paperconfig_{i}.yaml')
            with open(config_file, 'w') as f:
                yaml.dump(config, f)
            configs.append(config_file)

        # Create the list file
        list_file = os.path.join(temp_data_dir, 'paperconfig_files.txt')
        with open(list_file, 'w') as f:
            f.write('\n'.join(configs) + '\n')

        return list_file, configs

    @pytest.fixture
    def multi_adapter(self, two_config_list):
        """Create MultiOMICSAdapter and mock extracted_data on each inner adapter."""
        list_file, _ = two_config_list
        adapter = MultiOMICSAdapter(config_list_file=list_file)

        for i, inner in enumerate(adapter.adapters):
            inner.extracted_data = {
                'publication': {
                    'publication_id': f'test_pub_{i}',
                    'title': f'Test Publication {i}',
                    'doi': f'10.1234/test.{i}',
                },
            }
        return adapter

    def test_loads_all_configs(self, multi_adapter):
        """Verify all config files are loaded."""
        assert len(multi_adapter.adapters) == 2

    def test_get_nodes_combines_all(self, multi_adapter):
        """Verify get_nodes returns nodes from all configs."""
        nodes = multi_adapter.get_nodes()
        pub_nodes = [n for n in nodes if n[1] == 'publication']
        assert len(pub_nodes) == 2, f"Expected 2 publication nodes, got {len(pub_nodes)}"

        # OMICS adapter no longer creates organism nodes (CyanorakNcbi handles that)
        organism_nodes = [n for n in nodes if n[1] == 'organism']
        assert len(organism_nodes) == 0, f"Expected 0 organism nodes from OMICS, got {len(organism_nodes)}"

    def test_get_edges_combines_all(self, multi_adapter):
        """Verify get_edges returns expression edges from all configs."""
        edges = multi_adapter.get_edges()
        # 2 genes per config x 2 configs = 4 expression edges
        # Plus has_experiment edges (1 per experiment per publication)
        # Plus tests_coculture_with edges (1 per coculture experiment)
        expr_edges = [
            e for e in edges
            if e[3] == 'changes_expression_of'
        ]
        assert len(expr_edges) == 4, f"Expected 4 expression edges, got {len(expr_edges)}"

    def test_edges_have_distinct_sources(self, multi_adapter):
        """Verify expression edges come from different experiment sources."""
        edges = multi_adapter.get_edges()
        expr_edges = [
            e for e in edges
            if e[3] == 'changes_expression_of'
        ]
        source_ids = set(e[1] for e in expr_edges)
        assert len(source_ids) == 2, f"Expected 2 distinct experiment sources, got {source_ids}"

    def test_skips_blank_and_comment_lines(self, temp_data_dir):
        """Verify blank lines and comments in list file are skipped."""
        # Create a minimal config
        data_file = os.path.join(temp_data_dir, 'de.csv')
        pd.DataFrame({
            'Synonym': ['G1'], 'log2_fold_change': [1.0], 'adjusted_p_value': [0.01],
        }).to_csv(data_file, index=False)

        config = {
            'publication': {
                'papername': 'Test',
                'experiments': {
                    'test_exp': {
                        'name': 'Test',
                        'organism': 'Prochlorococcus MED4',
                        'treatment_condition': 'Treatment',
                        'control_condition': 'Control',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'coculture',
                        'treatment_organism': 'Org',
                        'treatment_taxid': 12345,
                    }
                },
                'supplementary_materials': {
                    's1': {
                        'type': 'csv', 'filename': data_file,
                        'statistical_analyses': [{
                            'id': 'test_list_adapter',
                            'experiment': 'test_exp',
                            'name_col': 'Synonym', 'logfc_col': 'log2_fold_change',
                            'adjusted_p_value_col': 'adjusted_p_value',
                        }]
                    }
                }
            }
        }
        config_file = os.path.join(temp_data_dir, 'pc.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        list_file = os.path.join(temp_data_dir, 'list.txt')
        with open(list_file, 'w') as f:
            f.write(f'# This is a comment\n\n{config_file}\n\n# Another comment\n')

        adapter = MultiOMICSAdapter(config_list_file=list_file)
        assert len(adapter.adapters) == 1

    def test_download_data_called_on_all(self, two_config_list):
        """Verify download_data delegates to all adapters."""
        list_file, _ = two_config_list
        adapter = MultiOMICSAdapter(config_list_file=list_file)
        adapter.download_data(cache=False)
        # After download_data, each adapter should have extracted_data attribute
        for inner in adapter.adapters:
            assert hasattr(inner, 'extracted_data') or hasattr(inner, 'pdf_extractor')


class TestGeneIdCleaning:
    """Test that gene IDs are cleaned (whitespace and asterisks stripped) before use."""

    def _make_adapter(self, temp_data_dir, data, analysis_extra=None):
        """Helper: create adapter from a DataFrame."""
        data_file = os.path.join(temp_data_dir, 'gene_id_clean.csv')
        data.to_csv(data_file, index=False)

        analysis = {
            'id': 'test_gene_id_clean',
            'experiment': 'gene_clean_exp',
            'name_col': 'Gene_ID',
            'logfc_col': 'log2FC',
            'adjusted_p_value_col': None,
        }
        if analysis_extra:
            analysis.update(analysis_extra)

        config = {
            'publication': {
                'papername': 'Test Gene ID Clean 2024',
                'experiments': {
                    'gene_clean_exp': {
                        'name': 'Test Gene ID Clean',
                        'organism': 'Prochlorococcus MED4',
                        'treatment_condition': 'Treatment',
                        'control_condition': 'Control',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'coculture',
                        'treatment_organism': 'Alteromonas macleodii',
                        'treatment_taxid': 28108,
                        'treatment_assembly_accession': 'GCF_001077695.1',
                    }
                },
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [analysis],
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_gene_clean.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'test_gene_clean',
                'title': 'Test',
                'doi': '10.1234/gene_clean',
            },
        }
        return adapter

    def test_gene_id_whitespace_stripped(self, temp_data_dir):
        """Verify leading/trailing whitespace in gene IDs is stripped."""
        data = pd.DataFrame({
            'Gene_ID': ['  GENE_A  ', 'GENE_B', '  GENE_C'],
            'log2FC': [1.0, -1.0, 0.5],
        })
        adapter = self._make_adapter(temp_data_dir, data)
        edges = adapter.get_edges()
        target_ids = [e[2] for e in edges if e[3] == 'changes_expression_of']
        target_str = ' '.join(target_ids)
        assert 'GENE_A' in target_str, "GENE_A (with spaces) should be found after stripping"
        assert 'GENE_B' in target_str
        assert 'GENE_C' in target_str
        assert len(target_ids) == 3

    def test_gene_id_asterisk_stripped(self, temp_data_dir):
        """Verify asterisks in gene IDs are stripped."""
        data = pd.DataFrame({
            'Gene_ID': ['GENE_A*', '*GENE_B', 'GENE_C'],
            'log2FC': [1.0, -1.0, 0.5],
        })
        adapter = self._make_adapter(temp_data_dir, data)
        edges = adapter.get_edges()
        target_ids = [e[2] for e in edges if e[3] == 'changes_expression_of']
        target_str = ' '.join(target_ids)
        assert 'GENE_A' in target_str, "GENE_A* should be cleaned to GENE_A"
        assert 'GENE_B' in target_str, "*GENE_B should be cleaned to GENE_B"
        assert 'GENE_C' in target_str
        assert len(target_ids) == 3

    def test_gene_id_whitespace_and_asterisk_stripped(self, temp_data_dir):
        """Verify both whitespace and asterisks are stripped from gene IDs."""
        data = pd.DataFrame({
            'Gene_ID': ['  GENE_A *  ', ' * GENE_B', 'GENE_C'],
            'log2FC': [1.0, -1.0, 0.5],
        })
        adapter = self._make_adapter(temp_data_dir, data)
        edges = adapter.get_edges()
        expr_edges = [e for e in edges if e[3] == 'changes_expression_of']
        assert len(expr_edges) == 3, f"Expected 3 expression edges after stripping, got {len(expr_edges)}"

    def test_gene_id_only_asterisks_skipped(self, temp_data_dir):
        """Verify gene IDs that reduce to empty string after stripping are skipped."""
        data = pd.DataFrame({
            'Gene_ID': ['***', 'GENE_B', '  *  '],
            'log2FC': [1.0, -1.0, 0.5],
        })
        adapter = self._make_adapter(temp_data_dir, data)
        edges = adapter.get_edges()
        # '***' and '  *  ' reduce to '' → should be skipped
        expr_edges = [e for e in edges if e[3] == 'changes_expression_of']
        assert len(expr_edges) == 1, f"Expected 1 expression edge (only GENE_B), got {len(expr_edges)}"


class TestUnifiedEdgeLabel:
    """Test that all expression edges use the unified changes_expression_of label."""

    def test_all_expression_edges_use_changes_expression_of(self, temp_data_dir):
        """Both coculture and non-coculture experiments use changes_expression_of."""
        data_file = os.path.join(temp_data_dir, 'de.csv')
        pd.DataFrame({
            'Gene': ['PMM0001', 'PMM0002'],
            'log2FC': [1.5, -2.0],
            'padj': [0.01, 0.05],
        }).to_csv(data_file, index=False)

        config = {
            'publication': {
                'papername': 'Test Unified Label 2024',
                'experiments': {
                    'light_exp': {
                        'name': 'High light stress',
                        'organism': 'Prochlorococcus MED4',
                        'treatment_condition': 'High light',
                        'control_condition': 'Low light',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'light_stress',
                    }
                },
                'supplementary_materials': {
                    'table1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [{
                            'id': 'unified_label_test',
                            'experiment': 'light_exp',
                            'name_col': 'Gene',
                            'logfc_col': 'log2FC',
                            'adjusted_p_value_col': 'padj',
                        }],
                    }
                },
            }
        }
        config_file = os.path.join(temp_data_dir, 'pc_unified_label.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {'publication_id': 'unified_pub', 'doi': '10.1/ul'},
        }
        edges = adapter.get_edges()
        expr_edges = [e for e in edges if e[3] == 'changes_expression_of']
        assert len(expr_edges) == 2
        for edge in expr_edges:
            assert edge[3] == 'changes_expression_of'

    def test_no_old_edge_labels(self, adapter_with_mock_extracted_data):
        """Verify old edge labels are never emitted."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()
        old_labels = [e for e in edges if e[3] in (
            'affects_expression_of', 'condition_changes_expression_of',
            'coculture_changes_expression_of', 'published_expression_data_about'
        )]
        assert len(old_labels) == 0,             f"Old edge labels should never be emitted, got {[e[3] for e in old_labels]}"


class TestExperimentNodeProperties:
    """Test experiment node properties (organism_name, omics_type, medium, temperature, etc.)."""

    def test_experiment_has_organism_name(self, adapter_with_mock_extracted_data):
        """organism_name is set on experiment nodes."""
        nodes = adapter_with_mock_extracted_data.get_nodes()
        exp_nodes = [n for n in nodes if n[1] == 'experiment']
        assert len(exp_nodes) == 1
        assert exp_nodes[0][2]['organism_name'] == 'Prochlorococcus MED4'

    def test_experiment_has_medium_and_temperature(self, adapter_with_mock_extracted_data):
        """medium and temperature are set on experiment nodes."""
        nodes = adapter_with_mock_extracted_data.get_nodes()
        exp_nodes = [n for n in nodes if n[1] == 'experiment']
        assert len(exp_nodes) == 1
        props = exp_nodes[0][2]
        assert props['medium'] == 'Pro99'
        assert props['temperature'] == '24C'

    def test_experiment_has_coculture_partner(self, adapter_with_mock_extracted_data):
        """coculture_partner is set from treatment_organism."""
        nodes = adapter_with_mock_extracted_data.get_nodes()
        exp_nodes = [n for n in nodes if n[1] == 'experiment']
        assert len(exp_nodes) == 1
        assert exp_nodes[0][2]['coculture_partner'] == 'Alteromonas macleodii'

    def test_experiment_has_table_scope(self, adapter_with_mock_extracted_data):
        """table_scope is set on experiment nodes."""
        nodes = adapter_with_mock_extracted_data.get_nodes()
        exp_nodes = [n for n in nodes if n[1] == 'experiment']
        assert len(exp_nodes) == 1
        assert exp_nodes[0][2]['table_scope'] == 'all_detected_genes'

    def test_experiment_table_scope_defaults_to_empty(self, temp_data_dir, sample_de_data):
        """table_scope defaults to empty string when not set in paperconfig."""
        data_file = os.path.join(temp_data_dir, 'de_genes2.csv')
        sample_de_data.to_csv(data_file, index=False)
        config = {
            'publication': {
                'papername': 'Test No Scope',
                'experiments': {
                    'test_exp': {
                        'name': 'Test',
                        'organism': 'Prochlorococcus MED4',
                        'treatment_condition': 'Treatment',
                        'control_condition': 'Control',
                        'omics_type': 'RNASEQ',
                        'test_type': 'DESeq2',
                        'treatment_type': 'nitrogen',
                    }
                },
                'supplementary_materials': {
                    'supp1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [{
                            'id': 'a1',
                            'experiment': 'test_exp',
                            'name_col': 'Synonym',
                            'logfc_col': 'log2_fold_change',
                            'adjusted_p_value_col': 'adjusted_p_value',
                        }]
                    }
                }
            }
        }
        import yaml
        config_file = os.path.join(temp_data_dir, 'paperconfig_no_scope.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)
        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {'title': 'T', 'authors': [], 'journal': '', 'year': 2024, 'abstract': '', 'description': ''}
        nodes = adapter.get_nodes()
        exp_nodes = [n for n in nodes if n[1] == 'experiment']
        assert exp_nodes[0][2]['table_scope'] == ''


class TestHasExperimentEdges:
    """Test has_experiment edges from Publication -> Experiment."""

    def test_has_experiment_edges_emitted(self, adapter_with_mock_extracted_data):
        """Verify has_experiment edges are emitted from Publication to Experiment."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()
        has_exp_edges = [e for e in edges if e[3] == 'has_experiment']

        assert len(has_exp_edges) == 1,             f"Expected 1 has_experiment edge, got {len(has_exp_edges)}"

        edge = has_exp_edges[0]
        _, source_id, target_id, label, props = edge
        assert 'doi' in source_id.lower() or '10.1234' in source_id
        assert 'coculture_alteromonas_med4_rnaseq' in target_id


class TestTestsCocultureWithEdges:
    """Test tests_coculture_with edges from Experiment -> OrganismTaxon."""

    def test_coculture_edge_emitted(self, adapter_with_mock_extracted_data):
        """Verify tests_coculture_with edge is emitted for coculture experiments."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()
        coculture_edges = [e for e in edges if e[3] == 'tests_coculture_with']

        assert len(coculture_edges) == 1,             f"Expected 1 tests_coculture_with edge, got {len(coculture_edges)}"

        edge = coculture_edges[0]
        _, source_id, target_id, label, props = edge
        assert 'coculture_alteromonas_med4_rnaseq' in source_id
        assert 'insdc.gcf' in target_id.lower() or 'GCF_001077695.1' in target_id


if __name__ == '__main__':
    pytest.main([__file__, '-v'])