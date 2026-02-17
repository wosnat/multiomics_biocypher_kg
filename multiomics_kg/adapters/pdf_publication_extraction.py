"""
PDF Publication Node Adapter

Extracts publication metadata from PDF files using LangChain LLM API calls.
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


import dotenv
dotenv.load_dotenv()  # Load environment variables from .env file

try:
    from pypdf import PdfReader
except Exception:
    PdfReader = None

from langchain.chat_models import init_chat_model
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import JsonOutputParser

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")




class PDFPublicationExtractor:
    """
    Adapter for extracting publication metadata from PDF files
    using LangChain with any supported LLM provider (OpenAI, Anthropic, etc.).
    Results are cached to avoid redundant API calls.
    """

    def __init__(
        self,
        #model_name: str = "claude-haiku-4-5",
        model_name: str = "gpt-5-nano",
        temperature: float = 0.0,
        cache_file: str | Path = None,
    ):
        """
        Initialize the PDF Publication Adapter.

        Args:
            model_name: LangChain model string (e.g., "openai:gpt-4-turbo", "anthropic:claude-3-5-sonnet-20241022")
            temperature: LLM temperature setting
            cache_file: Path to cache JSON file for extracted metadata

        Example:
            from langchain_openai import ChatOpenAI
            from multiomics_kg.adapters.pdf_publication_adapter import PDFPublicationAdapter

            llm = ChatOpenAI(model="gpt-4-turbo", temperature=0)
            extractor = PDFExtractor(
                llm=llm,
                cache_file="pdf_extraction_cache.json"
            )
        """

        self.model_name = model_name
        # example models:
        # o3_mini = init_chat_model("openai:o3-mini", temperature=0)
        # claude_sonnet = init_chat_model("anthropic:claude-sonnet-4-5-20250929", temperature=0)
        # gemini_2-5_flash = init_chat_model("google_vertexai:gemini-2.5-flash", temperature=0)
        self.llm = init_chat_model(model=model_name, temperature=temperature)
        self.cache_file = Path(cache_file) if cache_file else Path("cache/pdf_extraction_cache.json")
        self.cache = self._load_cache()
        self.id_counter = dict()  # Track IDs for uniqueness: {base_id: count}
        self._setup_extraction_chain()

    def _setup_extraction_chain(self):
        """Setup LangChain extraction chain with JSON parser."""

        self.extraction_prompt = ChatPromptTemplate.from_messages([
            ("system", "You are an expert scientific paper analyzer. Extract structured metadata from academic papers. Return ONLY valid JSON, no additional text. CRITICAL: Always extract DOI if present - it is essential for unique identification."),
            ("human", """{text}

Extract the following metadata in JSON format. CRITICAL: Always extract and include DOI if available:
{{
  "publication": {{
    "title": "Paper title",
    "authors": ["Author 1", "Author 2"],
    "journal": "Journal name",
    "publication_year": "YYYY",
    "doi": "DOI if available (e.g., 10.1038/nature12345)",
    "pmid": "PubMed ID if available",
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

    def _extract_pdf_text(self, pdf_path: Path, max_pages: int = 5) -> str:
        """
        Extract text from PDF file, stopping at the References section.

        Args:
            pdf_path: Path to PDF file
            max_pages: Maximum pages to extract (to avoid huge PDFs)

        Returns:
            Extracted text from PDF, excluding References and beyond
        """
        if PdfReader is None:
            logger.warning("No PDF reader installed. Install with: pip install pypdf")
            return ""

        try:
            text = ""
            with open(pdf_path, "rb") as f:
                reader = PdfReader(f)
                # pypdf/PyPDF2 provide pages as a sequence-like object
                num_pages = min(len(reader.pages), max_pages)

                for page_num in range(num_pages):
                    page = reader.pages[page_num]
                    # pypdf uses `extract_text()`; older PyPDF2 variants may differ
                    page_text = getattr(page, 'extract_text', lambda: '')()
                    
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

            return text[:10000]  # Limit to first 10000 chars to avoid token limits
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
        cache_key = str(pdf_path)  # Use original string before Path normalization
        pdf_path = Path(pdf_path)

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
        #pdf_text = self._clean_pdf_text(pdf_text)
        # print(f'After cleaning:')
        # print(f'pdf_text length: {len(pdf_text)}')
        # print(f'pdf_text sample: {pdf_text[:500]}') 

        #try:
        # Use LangChain chain to extract metadata
        extracted = self.chain.invoke({"text": pdf_text})

        # Generate and store publication ID based on DOI
        if "publication" in extracted:
            pub = extracted["publication"]
            # Prefer DOI for ID (most unique), then pmid, then generate from title
            pub_doi = pub.get("doi", "") or ""
            pub_doi = pub_doi.strip()
            pub_id = pub_doi or pub.get("pmid", "") or f"pub_{pub.get('title', 'unknown')[:30]}"
            pub_id = self._generate_unique_id(pub_id, "pub")
            extracted["publication"]["publication_id"] = pub_id

        # Cache the result
        self.cache[cache_key] = extracted
        self._save_cache()

        logger.info(f"Successfully extracted metadata from {pdf_path.name}")
        return extracted

        # except Exception as e:
        #     logger.error(f"LLM extraction failed for {pdf_path}: {e}")
        #     return {}


 



# Example usage
if __name__ == "__main__":
    import sys

    import dotenv
    dotenv.load_dotenv()  # Load environment variables from .env file
    adapter = PDFPublicationExtractor(
        cache_file="pdf_extraction_cache.json",
    )


    pdf_fpath = "data/Prochlorococcus/papers_and_supp/Aharonovich 2016/41396_2016_article_bfismej201670.pdf"
    pdf_fpath = "data/Prochlorococcus/papers_and_supp/Anjur 2025/2025.08.05.668435v1.full.pdf"
    pdf_text = adapter.extract_from_pdf(pdf_fpath)
    
    print(json.dumps(pdf_text, indent=2))
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
