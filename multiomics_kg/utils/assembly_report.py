"""NCBI assembly-quality report for a genome strain.

Step 0 (``download_genome_data._ncbi_fetch_assembly_reports``) writes one
``assembly_report.json`` per strain directory, shaped
``{"fetched_at": "YYYY-MM-DD", "report": <NCBI Datasets v2 dataset_report>}``.
"""
from __future__ import annotations

ASSEMBLY_REPORT_FILENAME = "assembly_report.json"
