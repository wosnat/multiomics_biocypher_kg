"""
Tests for the OMICS adapter organism to gene expression association edges.

Tests verify that the adapter correctly creates:
1. Publication nodes (with study metadata merged in)
2. Organism nodes from treatment_organism and treatment_taxid
3. Environmental condition nodes from environmental_conditions section
4. affects_expression_of edges from organism → gene
5. Correct edge properties (log2FC, p-value, direction, experimental_context, etc.)
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

    # Create config with environmental_conditions section
    config = {
        'publication': {
            'papername': 'Test Publication 2024',
            'environmental_conditions': {
                'test_growth_conditions': {
                    'condition_type': 'growth_medium',
                    'name': 'Test growth conditions',
                    'light_condition': 'continuous_light',
                    'temperature': '24C',
                }
            },
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
                            'control_condition': 'Axenic',
                            'treatment_condition': 'Coculture with Alteromonas',
                            'experimental_context': 'in test growth medium under continuous light',
                            'timepoint': '24h',
                            'organism': 'Prochlorococcus MED4',
                            'treatment_organism': 'Alteromonas macleodii',
                            'treatment_taxid': 28108,
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
        assert properties.get('organism') == ['Prochlorococcus MED4']


class TestOrganismNodeCreation:
    """Test that organism nodes are created correctly."""

    def test_organism_node_created_from_treatment_taxid(self, adapter_with_mock_extracted_data):
        """Verify organism node is created with correct taxid."""
        adapter = adapter_with_mock_extracted_data
        nodes = adapter.get_nodes()

        # Find organism nodes
        organism_nodes = [n for n in nodes if n[1] == 'organism']

        assert len(organism_nodes) == 1, "Expected exactly one organism node"

        org_node = organism_nodes[0]
        node_id, label, properties = org_node

        # Check node ID uses ncbitaxon prefix
        assert 'ncbitaxon:28108' in node_id.lower() or '28108' in node_id, \
            f"Expected organism node ID to contain taxid 28108, got {node_id}"

        # Check organism name property
        assert properties.get('organism_name') == 'Alteromonas macleodii', \
            f"Expected organism_name 'Alteromonas macleodii', got {properties.get('organism_name')}"

    def test_no_organism_node_without_treatment_taxid(self, temp_data_dir, sample_de_data):
        """Verify no organism node is created when treatment_taxid is missing."""
        # Save sample DE data
        data_file = os.path.join(temp_data_dir, 'de_genes.csv')
        sample_de_data.to_csv(data_file, index=False)

        # Config without treatment_taxid
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
                                # No treatment_organism or treatment_taxid
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

        assert len(organism_nodes) == 0, "No organism node should be created without treatment_taxid"


class TestEnvironmentalConditionNodeCreation:
    """Test that environmental condition nodes are created correctly."""

    def test_environmental_condition_node_created(self, adapter_with_mock_extracted_data):
        """Verify environmental condition node is created from config."""
        adapter = adapter_with_mock_extracted_data
        nodes = adapter.get_nodes()

        # Find environmental condition nodes
        env_nodes = [n for n in nodes if n[1] == 'environmental_condition']

        assert len(env_nodes) == 1, f"Expected exactly one environmental condition node, got {len(env_nodes)}"

        env_node = env_nodes[0]
        node_id, label, properties = env_node

        # Check that ID includes DOI prefix for uniqueness
        assert 'test_pub_2024' in node_id or '10.1234' in node_id, \
            f"Expected environmental condition node ID to include publication ID, got {node_id}"
        assert 'test_growth_conditions' in node_id, \
            f"Expected environmental condition node ID to include local ID, got {node_id}"

        # Check properties are copied from config
        assert properties.get('condition_type') == 'growth_medium', \
            f"Expected condition_type 'growth_medium', got {properties.get('condition_type')}"
        assert properties.get('light_condition') == 'continuous_light', \
            f"Expected light_condition 'continuous_light', got {properties.get('light_condition')}"
        assert properties.get('local_id') == 'test_growth_conditions', \
            f"Expected local_id 'test_growth_conditions', got {properties.get('local_id')}"

    def test_no_env_condition_node_without_section(self, temp_data_dir, sample_de_data):
        """Verify no environmental condition node is created when section is missing."""
        data_file = os.path.join(temp_data_dir, 'de_genes.csv')
        sample_de_data.to_csv(data_file, index=False)

        # Config without environmental_conditions section
        config = {
            'publication': {
                'papername': 'Test Publication 2024',
                # No environmental_conditions section
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
                                'treatment_organism': 'Alteromonas macleodii',
                                'treatment_taxid': 28108,
                                'name_col': 'Synonym',
                                'logfc_col': 'log2_fold_change',
                                'adjusted_p_value_col': 'adjusted_p_value',
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_no_env.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        adapter = OMICSAdapter(config_file=config_file)
        adapter.extracted_data = {
            'publication': {'publication_id': 'test', 'doi': '10.1234/test'},
        }

        nodes = adapter.get_nodes()
        env_nodes = [n for n in nodes if n[1] == 'environmental_condition']

        assert len(env_nodes) == 0, "No environmental condition node should be created without section"


class TestOrganismToGeneEdges:
    """Test that organism → gene expression edges are created correctly."""

    def test_affects_expression_of_edge_created(self, adapter_with_mock_extracted_data):
        """Verify affects_expression_of edges are created."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        # Find affects_expression_of edges
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        assert len(expression_edges) == 4, f"Expected 4 expression edges (one per gene), got {len(expression_edges)}"

    def test_edge_source_is_organism(self, adapter_with_mock_extracted_data):
        """Verify edge source is the treatment organism."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            # Source should be the organism (ncbitaxon:28108)
            assert 'ncbitaxon' in source_id.lower() or '28108' in source_id, \
                f"Expected source to be organism taxid, got {source_id}"

    def test_edge_target_is_gene(self, adapter_with_mock_extracted_data):
        """Verify edge target is a gene."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            # Target should be a gene (ncbigene:PMM000X)
            assert 'ncbigene' in target_id.lower() or 'PMM' in target_id, \
                f"Expected target to be gene ID, got {target_id}"

    def test_edge_has_log2_fold_change(self, adapter_with_mock_extracted_data):
        """Verify edges have log2_fold_change property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

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

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'adjusted_p_value' in properties, \
                f"Expected adjusted_p_value in edge properties, got {properties.keys()}"

    def test_edge_has_expression_direction(self, adapter_with_mock_extracted_data):
        """Verify edges have expression_direction property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

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

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            fc = properties['log2_fold_change']
            direction = properties['expression_direction']

            if fc > 0:
                assert direction == 'up', f"Positive fold change ({fc}) should have direction 'up', got '{direction}'"
            else:
                assert direction == 'down', f"Negative fold change ({fc}) should have direction 'down', got '{direction}'"

    def test_edge_has_control_condition(self, adapter_with_mock_extracted_data):
        """Verify edges have control_condition property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'control_condition' in properties, \
                f"Expected control_condition in edge properties"
            assert properties['control_condition'] == 'Axenic'

    def test_edge_has_experimental_context(self, adapter_with_mock_extracted_data):
        """Verify edges have experimental_context property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'experimental_context' in properties, \
                f"Expected experimental_context in edge properties"
            assert 'continuous light' in properties['experimental_context']

    def test_edge_has_time_point(self, adapter_with_mock_extracted_data):
        """Verify edges have time_point property when specified in config."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'time_point' in properties, \
                f"Expected time_point in edge properties"
            assert properties['time_point'] == '24h', \
                f"Expected time_point '24h', got {properties['time_point']}"

    def test_edge_has_publications(self, adapter_with_mock_extracted_data):
        """Verify edges have publications property."""
        adapter = adapter_with_mock_extracted_data
        edges = adapter.get_edges()

        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'publications' in properties, \
                f"Expected publications in edge properties"
            assert isinstance(properties['publications'], list), \
                f"Expected publications to be a list"


class TestEnvConditionToGeneEdges:
    """Test that environmental condition → gene expression edges are created correctly."""

    @pytest.fixture
    def env_condition_config(self, temp_data_dir, sample_de_data):
        """Create config with environmental_treatment_condition_id (no treatment_organism)."""
        data_file = os.path.join(temp_data_dir, 'de_genes.csv')
        sample_de_data.to_csv(data_file, index=False)

        config = {
            'publication': {
                'papername': 'Test Env Condition Publication 2024',
                'environmental_conditions': {
                    'test_gas_shock': {
                        'condition_type': 'gas_shock',
                        'name': 'Air (0.036% CO2 / 21% O2)',
                        'co2_level': '0.036%',
                        'oxygen_level': '21%',
                        'description': 'Standard air conditions',
                    }
                },
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [
                            {
                                'type': 'MICROARRAY',
                                'name': 'Gene expression response to air',
                                'id': 'gas_shock_air',
                                'test_type': 'Affymetrix microarray',
                                'control_condition': 'Time 0 (pre-shock)',
                                'treatment_condition': '0.036% CO2 / 21% O2 (air)',
                                'experimental_context': 'Pro99 medium, constant light',
                                'organism': 'Prochlorococcus MED4',
                                'environmental_treatment_condition_id': 'test_gas_shock',
                                'name_col': 'Synonym',
                                'logfc_col': 'log2_fold_change',
                                'adjusted_p_value_col': 'adjusted_p_value',
                            }
                        ]
                    }
                }
            }
        }

        config_file = os.path.join(temp_data_dir, 'paperconfig_env.yaml')
        with open(config_file, 'w') as f:
            yaml.dump(config, f)

        return config_file

    @pytest.fixture
    def env_adapter(self, env_condition_config):
        """Create adapter with env condition config and mocked extracted_data."""
        adapter = OMICSAdapter(config_file=env_condition_config)
        adapter.extracted_data = {
            'publication': {
                'publication_id': 'test_env_pub_2024',
                'title': 'Test Env Publication',
                'doi': '10.1234/test_env.2024',
                'authors': ['Test Author'],
                'description': 'Test gas shock study',
                'study_type': 'Transcriptomics',
                'organism': ['Prochlorococcus MED4'],
            },
        }
        return adapter

    def test_env_condition_edges_created(self, env_adapter):
        """Verify affects_expression_of edges are created with env condition source."""
        edges = env_adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']
        assert len(expression_edges) == 4, f"Expected 4 expression edges, got {len(expression_edges)}"

    def test_edge_source_is_env_condition(self, env_adapter):
        """Verify edge source is the environmental condition node ID."""
        edges = env_adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'test_gas_shock' in source_id, \
                f"Expected source to contain env condition ID 'test_gas_shock', got {source_id}"
            assert 'test_env_pub_2024' in source_id, \
                f"Expected source to contain publication ID, got {source_id}"

    def test_edge_target_is_gene(self, env_adapter):
        """Verify edge target is a gene."""
        edges = env_adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'ncbigene' in target_id.lower() or 'PMM' in target_id, \
                f"Expected target to be gene ID, got {target_id}"

    def test_edge_properties(self, env_adapter):
        """Verify edges have correct properties."""
        edges = env_adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'log2_fold_change' in properties
            assert 'adjusted_p_value' in properties
            assert 'expression_direction' in properties
            assert properties['expression_direction'] in ['up', 'down']
            assert 'control_condition' in properties
            assert properties['control_condition'] == 'Time 0 (pre-shock)'

    def test_no_organism_edges_with_env_condition(self, env_adapter):
        """Verify no organism-sourced edges when only env condition is specified."""
        edges = env_adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, target_id, label, properties = edge
            assert 'ncbitaxon' not in source_id.lower(), \
                f"Expected no organism source when env condition is specified, got {source_id}"


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
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': data_file,
                        'statistical_analyses': [
                            {
                                'type': 'RNASEQ',
                                'name': 'Test DE Analysis',
                                'id': 'test_missing_fc',
                                'test_type': 'DESeq2',
                                'control_condition': 'Control',
                                'treatment_condition': 'Treatment',
                                'organism': 'Prochlorococcus MED4',
                                'treatment_organism': 'Alteromonas macleodii',
                                'treatment_taxid': 28108,
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
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']
        # Only PMM0001 (2.5) and PMM0005 (-1.0) have valid FC values
        assert len(expression_edges) == 2, \
            f"Expected 2 edges (skipping empty, NaN, NA), got {len(expression_edges)}"

    def test_skipped_rows_do_not_appear_in_edges(self, adapter_with_missing_fc):
        """Verify that genes with missing FC are not in the edge targets."""
        edges = adapter_with_missing_fc.get_edges()
        target_ids = [e[2] for e in edges if e[3] == 'affects_expression_of']
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
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': csv_with_multi_row_header,
                        'statistical_analyses': [
                            {
                                'type': 'RNASEQ',
                                'name': 'Test DE with skip_rows',
                                'id': 'test_skip_rows',
                                'test_type': 'DESeq2',
                                'control_condition': 'Axenic',
                                'treatment_condition': 'Coculture',
                                'organism': 'Prochlorococcus NATL2A',
                                'treatment_organism': 'Alteromonas macleodii',
                                'treatment_taxid': 28108,
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
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']
        assert len(expression_edges) == 3, \
            f"Expected 3 edges (GENE_A, GENE_B, GENE_C), got {len(expression_edges)}"

    def test_skip_rows_correct_gene_ids(self, adapter_skip_rows):
        """Verify the correct gene IDs are parsed after skipping rows."""
        edges = adapter_skip_rows.get_edges()
        target_ids = [e[2] for e in edges if e[3] == 'affects_expression_of']
        target_str = ' '.join(target_ids)
        assert 'GENE_A' in target_str
        assert 'GENE_B' in target_str
        assert 'GENE_C' in target_str

    def test_skip_rows_correct_fold_change(self, adapter_skip_rows):
        """Verify fold change values are correctly parsed after skipping rows."""
        edges = adapter_skip_rows.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']
        fc_values = sorted([e[4]['log2_fold_change'] for e in expression_edges])
        assert fc_values == pytest.approx(sorted([-0.37, -0.212, 0.425]))

    def test_no_skip_rows_fails_gracefully(self, temp_data_dir, csv_with_multi_row_header):
        """Verify that without skip_rows, multi-row header CSV fails to find columns."""
        config = {
            'publication': {
                'papername': 'Test No Skip 2024',
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': csv_with_multi_row_header,
                        'statistical_analyses': [
                            {
                                'type': 'RNASEQ',
                                'id': 'test_no_skip_rows',
                                'treatment_organism': 'Alteromonas macleodii',
                                'treatment_taxid': 28108,
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
        # Without skip_rows, the column Gene_ID won't be found → no edges
        assert len(edges) == 0


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
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': csv_with_asterisks,
                        'statistical_analyses': [
                            {
                                'type': 'RNASEQ',
                                'name': 'Test DE with asterisk pvalue',
                                'id': 'test_asterisk',
                                'test_type': 'DESeq2',
                                'control_condition': 'Axenic',
                                'treatment_condition': 'Coculture',
                                'organism': 'Prochlorococcus NATL2A',
                                'treatment_organism': 'Alteromonas macleodii',
                                'treatment_taxid': 28108,
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
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']
        assert len(expression_edges) == 4, \
            f"Expected 4 edges, got {len(expression_edges)}"

    def test_fold_change_values_parsed_correctly(self, adapter_asterisk):
        """Verify asterisks are stripped and fold change is parsed as float."""
        edges = adapter_asterisk.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']
        fc_values = sorted([e[4]['log2_fold_change'] for e in expression_edges])
        assert fc_values == pytest.approx(sorted([-0.212, 0.425, -0.37, 0.15]))

    def test_significant_rows_get_low_pvalue(self, adapter_asterisk):
        """Verify asterisk-marked rows get adjusted_p_value < 0.5."""
        edges = adapter_asterisk.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, _, target_id, _, properties = edge
            if 'GENE_A' in target_id or 'GENE_C' in target_id:
                assert 'adjusted_p_value' in properties
                assert properties['adjusted_p_value'] < 0.5, \
                    f"Significant gene should have low p-value, got {properties['adjusted_p_value']}"

    def test_non_significant_rows_get_high_pvalue(self, adapter_asterisk):
        """Verify non-asterisk rows get adjusted_p_value = 1.0."""
        edges = adapter_asterisk.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, _, target_id, _, properties = edge
            if 'GENE_B' in target_id or 'GENE_D' in target_id:
                assert 'adjusted_p_value' in properties
                assert properties['adjusted_p_value'] == 1.0, \
                    f"Non-significant gene should have p-value 1.0, got {properties['adjusted_p_value']}"

    def test_expression_direction_correct(self, adapter_asterisk):
        """Verify expression direction is correct after asterisk stripping."""
        edges = adapter_asterisk.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, _, target_id, _, properties = edge
            fc = properties['log2_fold_change']
            direction = properties['expression_direction']
            if fc > 0:
                assert direction == 'up'
            else:
                assert direction == 'down'

    def test_without_asterisk_flag_values_are_skipped(self, temp_data_dir, csv_with_asterisks):
        """Verify that without pvalue_asterisk_in_logfc, asterisk values cause rows to be skipped."""
        config = {
            'publication': {
                'papername': 'Test No Asterisk Flag 2024',
                'supplementary_materials': {
                    'supp_table_1': {
                        'type': 'csv',
                        'filename': csv_with_asterisks,
                        'statistical_analyses': [
                            {
                                'type': 'RNASEQ',
                                'id': 'test_no_asterisk',
                                'treatment_organism': 'Alteromonas macleodii',
                                'treatment_taxid': 28108,
                                'name_col': 'Gene_ID',
                                'logfc_col': 'log2FC',
                                'adjusted_p_value_col': None,
                                # No pvalue_asterisk_in_logfc
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
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']
        # Only GENE_B (0.425) and GENE_D (0.15) have no asterisk → parseable as float
        assert len(expression_edges) == 2, \
            f"Expected 2 edges (asterisk values should fail float parse), got {len(expression_edges)}"


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
                'environmental_conditions': {
                    'biller_coculture': {
                        'condition_type': 'coculture',
                        'name': 'Co-culture with Alteromonas',
                    }
                },
                'supplementary_materials': {
                    'supp_table_2': {
                        'type': 'csv',
                        'filename': biller_style_csv,
                        'statistical_analyses': [
                            {
                                'type': 'RNASEQ',
                                'name': 'DE coculture vs axenic',
                                'id': 'de_biller_style',
                                'test_type': 'DESeq2',
                                'control_condition': 'Axenic',
                                'treatment_condition': 'Coculture',
                                'organism': 'Prochlorococcus NATL2A',
                                'environmental_treatment_condition_id': 'biller_coculture',
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
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']
        assert len(expression_edges) == 4

    def test_combined_correct_fold_changes(self, adapter_combined):
        """Verify fold change values are correct with combined features."""
        edges = adapter_combined.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']
        fc_values = sorted([e[4]['log2_fold_change'] for e in expression_edges])
        assert fc_values == pytest.approx(sorted([-0.425, -0.37, -0.212, 0.048]))

    def test_combined_significance_markers(self, adapter_combined):
        """Verify significant vs non-significant p-values with combined features."""
        edges = adapter_combined.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        significant_count = sum(
            1 for e in expression_edges if e[4].get('adjusted_p_value', 1.0) < 0.5
        )
        non_significant_count = sum(
            1 for e in expression_edges if e[4].get('adjusted_p_value', 1.0) >= 0.5
        )
        assert significant_count == 3, f"Expected 3 significant edges, got {significant_count}"
        assert non_significant_count == 1, f"Expected 1 non-significant edge, got {non_significant_count}"

    def test_combined_edge_source_is_env_condition(self, adapter_combined):
        """Verify edge source is env condition (not organism) with combined features."""
        edges = adapter_combined.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        for edge in expression_edges:
            _, source_id, _, _, _ = edge
            assert 'biller_coculture' in source_id


class TestEdgeWithRealData:
    """Test with the real Aharonovich 2016 config if available."""

    @pytest.fixture
    def real_config_path(self):
        """Path to real config file."""
        return 'data/Prochlorococcus/papers_and_supp/Aharonovich 2016/paperconfig.yaml'

    def test_real_config_creates_organism_nodes(self, real_config_path):
        """Test with real Aharonovich 2016 config."""
        if not os.path.exists(real_config_path):
            pytest.skip("Real config file not found")

        adapter = OMICSAdapter(config_file=real_config_path)
        adapter.download_data(cache=True)

        nodes = adapter.get_nodes()
        organism_nodes = [n for n in nodes if n[1] == 'organism']

        # Should have at least one organism node for Alteromonas
        assert len(organism_nodes) >= 1, "Expected at least one organism node"

        # Check for Alteromonas taxid
        org_ids = [n[0] for n in organism_nodes]
        has_alteromonas = any('28108' in str(oid) for oid in org_ids)
        assert has_alteromonas, f"Expected Alteromonas (taxid 28108) in organism nodes, got {org_ids}"

    def test_real_config_creates_environmental_condition_nodes(self, real_config_path):
        """Test that real config creates environmental condition nodes."""
        if not os.path.exists(real_config_path):
            pytest.skip("Real config file not found")

        adapter = OMICSAdapter(config_file=real_config_path)
        adapter.download_data(cache=True)

        nodes = adapter.get_nodes()
        env_nodes = [n for n in nodes if n[1] == 'environmental_condition']

        # Should have at least one environmental condition node
        assert len(env_nodes) >= 1, "Expected at least one environmental condition node"

        # Check that node has expected properties
        env_node = env_nodes[0]
        node_id, label, properties = env_node
        assert 'condition_type' in properties or 'name' in properties, \
            f"Expected condition_type or name in env node properties, got {properties.keys()}"

    def test_real_config_creates_expression_edges(self, real_config_path):
        """Test that real config creates expression edges."""
        if not os.path.exists(real_config_path):
            pytest.skip("Real config file not found")

        adapter = OMICSAdapter(config_file=real_config_path, test_mode=True)
        adapter.download_data(cache=True)

        edges = adapter.get_edges()
        expression_edges = [e for e in edges if e[3] == 'affects_expression_of']

        # Should have expression edges
        assert len(expression_edges) > 0, "Expected expression edges from real config"

        # Verify structure of first edge
        first_edge = expression_edges[0]
        _, source_id, target_id, label, properties = first_edge

        assert label == 'affects_expression_of'
        assert 'ncbitaxon' in source_id.lower() or '28108' in source_id
        assert 'log2_fold_change' in properties
        assert 'expression_direction' in properties


class TestMultiOMICSAdapter:
    """Test that MultiOMICSAdapter correctly aggregates results from multiple configs."""

    @pytest.fixture
    def two_config_list(self, temp_data_dir):
        """Create two separate paperconfig.yaml files and a list file pointing to them."""
        configs = []
        for i, (paper, organism, taxid) in enumerate([
            ('Paper A 2024', 'Alteromonas macleodii', 28108),
            ('Paper B 2024', 'Synechococcus sp.', 32049),
        ]):
            # Create data file for this config
            data_file = os.path.join(temp_data_dir, f'de_genes_{i}.csv')
            pd.DataFrame({
                'Synonym': [f'GENE_{i}_1', f'GENE_{i}_2'],
                'log2_fold_change': [1.5, -2.0],
                'adjusted_p_value': [0.01, 0.05],
            }).to_csv(data_file, index=False)

            config = {
                'publication': {
                    'papername': paper,
                    'supplementary_materials': {
                        'supp_table_1': {
                            'type': 'csv',
                            'filename': data_file,
                            'statistical_analyses': [
                                {
                                    'type': 'RNASEQ',
                                    'id': f'test_multi_adapter_{organism.replace(" ", "_")}',
                                    'test_type': 'DESeq2',
                                    'control_condition': 'Control',
                                    'treatment_condition': 'Treatment',
                                    'organism': 'Prochlorococcus MED4',
                                    'treatment_organism': organism,
                                    'treatment_taxid': taxid,
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

        organism_nodes = [n for n in nodes if n[1] == 'organism']
        assert len(organism_nodes) == 2, f"Expected 2 organism nodes, got {len(organism_nodes)}"

    def test_get_edges_combines_all(self, multi_adapter):
        """Verify get_edges returns edges from all configs."""
        edges = multi_adapter.get_edges()
        # 2 genes per config × 2 configs = 4 edges
        assert len(edges) == 4, f"Expected 4 edges, got {len(edges)}"

    def test_edges_have_distinct_sources(self, multi_adapter):
        """Verify edges come from different organism sources."""
        edges = multi_adapter.get_edges()
        source_ids = set(e[1] for e in edges)
        assert len(source_ids) == 2, f"Expected 2 distinct sources, got {source_ids}"

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
                'supplementary_materials': {
                    's1': {
                        'type': 'csv', 'filename': data_file,
                        'statistical_analyses': [{
                            'type': 'RNASEQ', 'id': 'test_list_adapter',
                            'test_type': 'DESeq2',
                            'treatment_organism': 'Org', 'treatment_taxid': 12345,
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


if __name__ == '__main__':
    pytest.main([__file__, '-v'])