"""
PDF Publication and Study Node Adapter

Extracts publication and study metadata from PDF files using LangChain LLM API calls.
Supports multiple LLM providers (OpenAI, Anthropic, etc.).
Caches results in JSON to avoid redundant LLM calls.

Author: AI Assistant
Date: 2025-01-20
"""

import json
import logging
import os
from datetime import datetime
from enum import Enum, auto
from pathlib import Path
from typing import Optional, Any
import re

try:
    import PyPDF2
except ImportError:
    PyPDF2 = None

try:
    from langchain_core.language_models import BaseLLM
    from langchain_core.prompts import ChatPromptTemplate
    from langchain_core.output_parsers import JsonOutputParser
except ImportError:
    BaseLLM = None

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class PublicationNodeType(Enum):
    """Define types of nodes the adapter can provide."""
    PUBLICATION = auto()
    STUDY = auto()


class PublicationField(Enum):
    """Fields for publication nodes."""
    PUBMED_ID = "pubmed_id"
    TITLE = "title"
    AUTHORS = "authors"
    JOURNAL = "journal"
    PUBLICATION_DATE = "publication_date"
    DOI = "doi"


class StudyField(Enum):
    """Fields for study nodes."""
    STUDY_ID = "study_id"
    TITLE = "title"
    DESCRIPTION = "description"
    ABSTRACT = "abstract"
    STUDY_TYPE = "study_type"
    ORGANISM = "organism"


class PDFPublicationAdapter:
    """
    Adapter for extracting publication and study metadata from PDF files
    using LangChain with any supported LLM provider (OpenAI, Anthropic, etc.).
    Results are cached to avoid redundant API calls.
    """

    def __init__(
        self,
        llm: BaseLLM,
        pdf_directory: str | Path = None,
        cache_file: str | Path = None,
    ):
        """
        Initialize the PDF Publication Adapter.

        Args:
            llm: LangChain BaseLLM instance (e.g., ChatOpenAI, ChatAnthropic)
            pdf_directory: Directory containing PDF files
            cache_file: Path to cache JSON file for extracted metadata

        Example:
            from langchain_openai import ChatOpenAI
            from multiomics_kg.adapters.pdf_publication_adapter import PDFPublicationAdapter

            llm = ChatOpenAI(model="gpt-4-turbo", temperature=0)
            adapter = PDFPublicationAdapter(
                llm=llm,
                pdf_directory="data/Prochlorococcus/papers_and_supp",
                cache_file="pdf_extraction_cache.json"
            )
        """
        if BaseLLM is None:
            raise ImportError(
                "LangChain not installed. Install with: pip install langchain"
            )

        self.llm = llm
        self.pdf_directory = Path(pdf_directory) if pdf_directory else None
        self.cache_file = Path(cache_file) if cache_file else Path("pdf_extraction_cache.json")
        self.cache = self._load_cache()
        self.id_counter = {}  # Track IDs for uniqueness: {base_id: count}
        self._setup_extraction_chain()

    def _setup_extraction_chain(self):
        """Setup LangChain extraction chain with JSON parser."""
        if BaseLLM is None:
            return

        self.extraction_prompt = ChatPromptTemplate.from_messages([
            ("system", "You are an expert scientific paper analyzer. Extract structured metadata from academic papers. Return ONLY valid JSON, no additional text. CRITICAL: Always extract DOI if present - it is essential for unique identification."),
            ("human", """{text}

Extract the following metadata in JSON format. CRITICAL: Always extract and include DOI if available:
{{
  "publication": {{
    "title": "Paper title",
    "authors": ["Author 1", "Author 2"],
    "journal": "Journal name",
    "publication_date": "YYYY-MM-DD or year if date unavailable",
    "doi": "DOI if available (e.g., 10.1038/nature12345)",
    "pubmed_id": "PubMed ID if available"
  }},
  "study": {{
    "title": "Study title (may be same as paper title)",
    "description": "Brief description of the study design and methods",
    "abstract": "Abstract or summary",
    "study_type": "e.g., Transcriptomics, Proteomics, Metabolomics, RNA-seq, Microarray",
    "organism": ["Organism 1", "Organism 2"]
  }}
}}

Return ONLY the JSON object.""")
        ])

        self.output_parser = JsonOutputParser()
        self.chain = self.extraction_prompt | self.llm | self.output_parser

    def _load_cache(self) -> dict:
        """Load extraction cache from JSON file."""
        if self.cache_file.exists():
            try:
                with open(self.cache_file) as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load cache: {e}")
        return {}

    def _save_cache(self):
        """Save extraction cache to JSON file."""
        try:
            with open(self.cache_file, "w") as f:
                json.dump(self.cache, f, indent=2, default=str)
        except Exception as e:
            logger.error(f"Failed to save cache: {e}")

    def _generate_unique_id(self, base_id: str, id_type: str = "pub") -> str:
        """
        Generate unique ID with counter for duplicates.

        Args:
            base_id: Base identifier (e.g., DOI or title)
            id_type: Type of ID ("pub" or "study")

        Returns:
            Unique ID string
        """
        key = f"{id_type}_{base_id}"
        if key not in self.id_counter:
            self.id_counter[key] = 0
            return base_id
        else:
            self.id_counter[key] += 1
            return f"{base_id}_v{self.id_counter[key]}"

    def _extract_pdf_text(self, pdf_path: Path, max_pages: int = 15) -> str:
        """
        Extract text from PDF file, stopping at the References section.

        Args:
            pdf_path: Path to PDF file
            max_pages: Maximum pages to extract (to avoid huge PDFs)

        Returns:
            Extracted text from PDF, excluding References and beyond
        """
        if PyPDF2 is None:
            logger.warning("PyPDF2 not installed. Install with: pip install PyPDF2")
            return ""

        try:
            text = ""
            with open(pdf_path, "rb") as f:
                reader = PyPDF2.PdfReader(f)
                num_pages = min(len(reader.pages), max_pages)

                for page_num in range(num_pages):
                    page = reader.pages[page_num]
                    page_text = page.extract_text()
                    
                    # Check if this page contains the References section
                    # Common patterns: "References", "REFERENCES", "References Cited", "Bibliography"
                    references_patterns = [
                        r'\nReferences\n',
                        r'\nREFERENCES\n',
                        r'\nBibliography\n',
                        r'\nBIBLIOGRAPHY\n',
                        r'^References\n',
                        r'^REFERENCES\n',
                    ]
                    
                    # Find the earliest occurrence of any reference pattern
                    references_index = -1
                    for pattern in references_patterns:
                        match = re.search(pattern, page_text, re.IGNORECASE)
                        if match:
                            references_index = match.start()
                            break
                    
                    if references_index != -1:
                        # Stop at references section
                        text += page_text[:references_index]
                        logger.debug(f"Stopped at References section on page {page_num}")
                        break
                    else:
                        text += page_text

            return text[:1000000]  # Limit to first 1000000 chars to avoid token limits
        except Exception as e:
            logger.error(f"Failed to extract text from {pdf_path}: {e}")
            return ""

    def _clean_pdf_text(self, text: str) -> str:
        """
        Clean extracted PDF text for LLM processing.

        Removes:
        - Excessive whitespace and formatting artifacts
        - Page numbers and headers/footers
        - Copyright notices and watermarks
        - Affiliation superscript numbers
        - Special characters and corrupted text

        Args:
            text: Raw extracted PDF text

        Returns:
            Cleaned text ready for LLM processing
        """
        # 1. Remove excessive whitespace
        text = re.sub(r'\n{3,}', '\n\n', text)  # Multiple newlines → double
        text = re.sub(r'[ \t]{2,}', ' ', text)  # Multiple spaces/tabs → single
        
        # 2. Remove page numbers and footers (common patterns)
        text = re.sub(r'^\s*\d+\s*$', '', text, flags=re.MULTILINE)  # Line with only number
        text = re.sub(r'^\s*[A-Z][a-z\s]+\s+\|\s+\d+\s*$', '', text, flags=re.MULTILINE)  # "Journal Name | 123"
        
        # 3. Remove common watermarks
        text = re.sub(r'(DRAFT|CONFIDENTIAL|Preprint|PREPRINT)\s*', '', text, flags=re.IGNORECASE)
        
        # 4. Remove copyright notices
        text = re.sub(r'©.*?\d{4}.*?(?=\n)', '', text, flags=re.IGNORECASE)
        text = re.sub(r'Copyright.*?\d{4}.*?(?=\n)', '', text, flags=re.IGNORECASE)
        
        # 5. Clean up superscript numbers after author names (e.g., "Smith1,2" → "Smith")
        #text = re.sub(r'([A-Za-z])[,\s]*\d+(?=[,\s])', r'\1', text)
        #text = re.sub(r'([A-Za-z])\d+(?=\n)', r'\1', text)
        
        # 6. Remove special characters that are likely OCR artifacts
        text = re.sub(r'[^\x20-\x7E\n\r]', '', text)  # Keep printable ASCII + newlines
        
        # 7. Clean up repeated special characters
        text = re.sub(r'[-_*]{3,}', '', text)  # Remove lines of dashes/underscores/asterisks
        
        # 8. Remove isolated numbers at line starts (likely page artifacts)
        text = re.sub(r'^\s*\d{1,3}\s*\n', '\n', text, flags=re.MULTILINE)
        
        # 9. Normalize spacing around common punctuation
        text = re.sub(r'\s+([.,;:?!])', r'\1', text)
        
        # 10. Remove leading/trailing whitespace from each line
        lines = [line.rstrip() for line in text.split('\n')]
        text = '\n'.join(lines)
        
        return text.strip()

    def extract_from_pdf(self, pdf_path: str | Path) -> dict:
        """
        Extract publication and study metadata from a PDF file using LangChain LLM.

        Args:
            pdf_path: Path to PDF file

        Returns:
            Dictionary with publication and study metadata
        """
        pdf_path = Path(pdf_path)
        cache_key = str(pdf_path)

        # Check cache first
        if cache_key in self.cache:
            # print(f'found in cache {cache_key}')
            logger.debug(f"Returning cached extraction for {pdf_path.name}")
            return self.cache[cache_key]

        # print(f'extracting {pdf_path.name}')
        logger.info(f"Extracting metadata from {pdf_path.name}...")

        # Extract text from PDF
        pdf_text = self._extract_pdf_text(pdf_path)
        # print(f'pdf_text length: {len(pdf_text)}')
        # print(f'pdf_text sample: {pdf_text[:500]}') 
        if not pdf_text:
            logger.warning(f"No text extracted from {pdf_path}")
            return {}

        # Clean the extracted text
        pdf_text = self._clean_pdf_text(pdf_text)
        # print(f'After cleaning:')
        # print(f'pdf_text length: {len(pdf_text)}')
        # print(f'pdf_text sample: {pdf_text[:500]}') 

        try:
            # Use LangChain chain to extract metadata
            extracted = self.chain.invoke({"text": pdf_text})

            # Generate and store IDs based on DOI
            if "publication" in extracted:
                pub = extracted["publication"]
                # Prefer DOI for ID (most unique), then pubmed_id, then generate from title
                pub_doi = pub.get("doi", "").strip()
                pub_id = pub_doi or pub.get("pubmed_id") or f"pub_{pub.get('title', 'unknown')[:30]}"
                pub_id = self._generate_unique_id(pub_id, "pub")
                extracted["publication"]["publication_id"] = pub_id

            if "study" in extracted:
                study = extracted["study"]
                # Prefer DOI from publication if available, otherwise generate from study title
                study_doi = None
                if "publication" in extracted:
                    study_doi = extracted["publication"].get("doi", "").strip()
                study_id = study_doi or f"study_{Path(pdf_path).stem}"
                study_id = self._generate_unique_id(study_id, "study")
                extracted["study"]["study_id"] = study_id

            # Cache the result
            self.cache[cache_key] = extracted
            self._save_cache()

            logger.info(f"Successfully extracted metadata from {pdf_path.name}")
            return extracted

        except Exception as e:
            logger.error(f"LLM extraction failed for {pdf_path}: {e}")
            return {}

    def extract_from_directory(self, directory: str | Path = None) -> list[dict]:
        """
        Extract metadata from all PDF files in a directory.

        Args:
            directory: Directory to search (uses self.pdf_directory if not provided)

        Returns:
            List of extracted metadata dictionaries
        """
        directory = Path(directory or self.pdf_directory)
        if not directory.exists():
            logger.error(f"Directory not found: {directory}")
            return []

        pdf_files = list(directory.glob("**/*.pdf"))
        logger.info(f"Found {len(pdf_files)} PDF files in {directory}")

        results = []
        for pdf_path in pdf_files:
            extracted = self.extract_from_pdf(pdf_path)
            if extracted:
                extracted["pdf_path"] = str(pdf_path)
                results.append(extracted)

        return results

    def get_nodes(self, node_type: PublicationNodeType = None) -> list[tuple]:
        """
        Generate BioCypher nodes from cached extractions.

        Args:
            node_type: Type of node to generate (PUBLICATION or STUDY)

        Returns:
            List of tuples (id, type, properties) for BioCypher
        """
        nodes = []

        for pdf_path, data in self.cache.items():
            if "publication" in data and (node_type is None or node_type == PublicationNodeType.PUBLICATION):
                pub = data["publication"]
                # Use stored ID from cache (which uses DOI if available)
                pub_id = pub.get("publication_id") or pub.get("doi") or pub.get("pubmed_id") or f"pub_{pub.get('title', 'unknown')[:20]}"

                nodes.append(
                    (
                        f"publication:{pub_id}",
                        "publication",
                        {
                            "title": pub.get("title"),
                            "authors": pub.get("authors", []),
                            "journal": pub.get("journal"),
                            "publication_date": pub.get("publication_date"),
                            "doi": pub.get("doi"),
                            "pubmed_id": pub.get("pubmed_id"),
                        },
                    )
                )

            if "study" in data and (node_type is None or node_type == PublicationNodeType.STUDY):
                study = data["study"]
                # Use stored ID from cache (which uses DOI if available)
                study_id = study.get("study_id") or f"study_{Path(pdf_path).stem}"

                nodes.append(
                    (
                        f"study:{study_id}",
                        "study",
                        {
                            "title": study.get("title"),
                            "description": study.get("description"),
                            "abstract": study.get("abstract"),
                            "study_type": study.get("study_type"),
                            "organism": study.get("organism", []),
                        },
                    )
                )

        return nodes

    def get_edges(self) -> list[tuple]:
        """
        Generate study_published_in edges linking studies to publications.

        Returns:
            List of tuples (source_id, target_id, edge_type)
        """
        edges = []

        for pdf_path, data in self.cache.items():
            if "publication" in data and "study" in data:
                pub = data["publication"]
                study = data["study"]

                # Use stored IDs from cache (which use DOI if available)
                pub_id = pub.get("publication_id") or pub.get("doi") or pub.get("pubmed_id") or f"pub_{pub.get('title', 'unknown')[:20]}"
                study_id = study.get("study_id") or f"study_{Path(pdf_path).stem}"

                edges.append(
                    (
                        f"study:{study_id}",
                        f"publication:{pub_id}",
                        "study_published_in",
                    )
                )

        return edges


# Example usage
if __name__ == "__main__":
    import sys

    import dotenv
    dotenv.load_dotenv()  # Load environment variables from .env file
    from langchain_openai import ChatOpenAI
    llm = ChatOpenAI(model="gpt-4-turbo", temperature=0)
    adapter = PDFPublicationAdapter(
        llm=llm,
        pdf_directory="data/Prochlorococcus/papers_and_supp",
        cache_file="pdf_extraction_cache.json",
    )


    pdf_fpath = "data/Prochlorococcus/papers_and_supp/Aharonovich 2016/41396_2016_article_bfismej201670.pdf"
    pdf_fpath = "data/Prochlorococcus/papers_and_supp/Anjur 2025/2025.08.05.668435v1.full.pdf"
    adapter.extract_from_pdf(pdf_fpath)
    nodes = adapter.get_nodes()
    edges = adapter.get_edges()

    with open('nodes.json', 'w', encoding='utf-8') as f:
        json.dump(nodes, f, indent=2, default=str)
    with open('edges.json', 'w', encoding='utf-8') as f:
        json.dump(edges, f, indent=2, default=str)
    
    # pdf_text = adapter._extract_pdf_text(Path(pdf_fpath))
    # pdf_text = adapter._clean_pdf_text(pdf_text)

    # with open('tmp.txt', 'w', encoding='utf-8') as f:
    #     f.write(pdf_text)   
    # #print(pdf_text)


#     logging.basicConfig(level=logging.INFO)


#     # Example with OpenAI
#     print("=" * 80)
#     print("Example usage with OpenAI:")
#     print("=" * 80)
#     print("""
# from langchain_openai import ChatOpenAI
# from multiomics_kg.adapters.pdf_publication_adapter import PDFPublicationAdapter

# # Initialize with OpenAI
# llm = ChatOpenAI(model="gpt-4-turbo", temperature=0, api_key="your-key")
# adapter = PDFPublicationAdapter(
#     llm=llm,
#     pdf_directory="data/Prochlorococcus/papers_and_supp",
#     cache_file="pdf_extraction_cache.json",
# )

# # Extract from a single PDF
# result = adapter.extract_from_pdf("path/to/paper.pdf")
# print(result)

# # Extract from directory
# results = adapter.extract_from_directory()

# # Get nodes and edges
# nodes = adapter.get_nodes()
# edges = adapter.get_edges()
# """)

#     print("\n" + "=" * 80)
#     print("Example usage with Anthropic:")
#     print("=" * 80)
#     print("""
# from langchain_anthropic import ChatAnthropic
# from multiomics_kg.adapters.pdf_publication_adapter import PDFPublicationAdapter

# # Initialize with Anthropic
# llm = ChatAnthropic(model="claude-3-5-sonnet-20241022", api_key="your-key")
# adapter = PDFPublicationAdapter(
#     llm=llm,
#     pdf_directory="data/Prochlorococcus/papers_and_supp",
#     cache_file="pdf_extraction_cache.json",
# )
# """)

#     print("\n" + "=" * 80)
#     print("Supported LangChain providers:")
#     print("=" * 80)
#     print("""
# - OpenAI: from langchain_openai import ChatOpenAI
# - Anthropic: from langchain_anthropic import ChatAnthropic
# - Google: from langchain_google_genai import ChatGoogleGenerativeAI
# - Cohere: from langchain_cohere import ChatCohere
# - Azure: from langchain_openai import AzureChatOpenAI
# - Ollama (local): from langchain_ollama import ChatOllama
# - And many more via LangChain integrations
# """)
