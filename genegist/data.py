import json
import os
from pathlib import Path
from typing import Sequence, Optional
import xml.etree.ElementTree as ET

from Bio import Entrez
import pandas as pd


class GeneRIFS:
    def __init__(self):
        """
        GeneRIFs (Gene References Into Function) are short, descriptive summaries of a
        gene's function from the Gene database at NCBI.
        """

        self.df = self.get_generifs()

    def get_generifs(self) -> pd.DataFrame:
        """Get the GeneRIFs data."""

        url = "https://ftp.ncbi.nlm.nih.gov/gene/GeneRIF/generifs_basic.gz"
        cache = Path(__file__).parent / "data" / "generifs_basic.parquet"
        if cache.exists():
            df = pd.read_parquet(cache)
        else:
            df = pd.read_csv(url, compression="gzip", sep="\t", lineterminator="\n")

            # Ensure that the 'PubMed ID (PMID) list' column is treated as a string
            df["PubMed ID (PMID) list"] = df["PubMed ID (PMID) list"].astype(str)

            cache.parent.mkdir(exist_ok=True)
            df.to_parquet(cache)
        return df

    def get_texts_by_gene(self, gene_name: str) -> Sequence[str]:
        """Get all the texts for the given gene IDs."""

        if gene_name.isnumeric():
            gene_id = int(gene_name)
        else:
            gene_id = gene2id(gene_name)
        mask = (self.df["#Tax ID"] == 9606) & (self.df["Gene ID"] == gene_id)
        return self.df[mask]["GeneRIF text"].tolist()

    def get_texts_by_gene_set(self, gene_set: str, add_abstracts: bool = False) -> dict:
        """Get all the texts for the given gene set."""

        gene_symbols = geneset2symbols(gene_set)
        texts = {}
        for gene_symbol in gene_symbols:
            texts[gene_symbol] = self.get_texts_by_gene(gene_symbol, add_abstracts)
        return texts

    def get_texts_by_gene_set_list(
        self, gene_set: list, add_abstracts: bool = False
    ) -> dict:
        """Get all the texts for the given gene set."""

        texts = {}
        for gene_symbol in gene_set:
            texts[gene_symbol] = self.get_texts_by_gene(gene_symbol, add_abstracts)
        return texts


def set_ncbi_email() -> None:
    """Set the NCBI email address."""

    if os.environ.get("NCBI_EMAIL"):
        Entrez.email = os.environ["NCBI_EMAIL"]
    else:
        raise ValueError("Please set the NCBI_EMAIL environment variable.")


def gene2id(gene_name: str) -> int:
    """Get the NCBI gene ID for the given gene name."""

    set_ncbi_email()

    handle = Entrez.esearch(
        db="gene", term=f"{gene_name}[Gene Name] AND Homo sapiens[Organism]"
    )
    record = Entrez.read(handle)
    handle.close()

    # Check if there are any results
    if not record["IdList"]:
        return None

    gene_id = record["IdList"][0]

    return int(gene_id)


def id2gene(gene_id: int) -> str:
    """Get the gene name for the given NCBI gene ID."""

    set_ncbi_email()

    handle = Entrez.esummary(db="gene", id=str(gene_id))
    record = Entrez.read(handle)
    handle.close()

    if not record:
        return None

    # Extracting the gene name from the record
    if (
        "DocumentSummarySet" in record
        and "DocumentSummary" in record["DocumentSummarySet"]
    ):
        summaries = record["DocumentSummarySet"]["DocumentSummary"]
        if summaries:
            gene_name = summaries[0]["Name"]
            return gene_name

    return None


def get_gene_abstracts(gene_name: str, max_results: int = 100) -> str:
    """Get NCBI published abstractions for a given gene."""

    set_ncbi_email()

    search_term = f"{gene_name}[Gene Name] AND Homo sapiens[Organism]"
    search_handle = Entrez.esearch(db="pubmed", term=search_term, retmax=max_results)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    id_list = search_results["IdList"]

    if id_list:
        fetch_handle = Entrez.efetch(
            db="pubmed", id=id_list, rettype="abstract", retmode="text"
        )
        abstracts = fetch_handle.read()
        fetch_handle.close()
        return abstracts
    else:
        return "No results found."


def get_abstract(pmid: str) -> str:
    """Get NCBI published abstract for a given PubMed ID."""

    set_ncbi_email()

    if pmid:
        fetch_handle = Entrez.efetch(
            db="pubmed", id=pmid, rettype="abstract", retmode="text"
        )
        abstract = fetch_handle.read()
        fetch_handle.close()
        return abstract
    else:
        return "No PMID provided."


def get_article(pmid: str) -> Optional[str]:
    """Get article for a given PMID, assuming it is open access."""

    set_ncbi_email()

    def parse_article_xml(xml_data: str) -> Optional[str]:
        try:
            root = ET.fromstring(xml_data)
            article_text_segments = []
            tags_of_interest = ["title", "p", "sec"]

            def extract_text(element):
                for el in element:
                    if el.tag in tags_of_interest:
                        if el.text:
                            article_text_segments.append(el.text.strip())
                    extract_text(el)

            extract_text(root)
            return "\n\n".join(article_text_segments)
        except Exception as e:
            return f"An error occurred while parsing XML: {e}"

    try:
        link_handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        link_results = Entrez.read(link_handle)
        link_handle.close()

        if not link_results[0]["LinkSetDb"]:
            return "No PMC article found for this PMID."

        pmc_id = link_results[0]["LinkSetDb"][0]["Link"][0]["Id"]

        fetch_handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml")
        article_xml = fetch_handle.read().strip()
        fetch_handle.close()

        return parse_article_xml(article_xml)
    except Exception as e:
        return f"An error occurred while fetching the article: {e}"


def geneset2symbols(geneset: str) -> Sequence[str]:
    """Get the gene symbols for the given geneset."""

    hallmarks_path = Path(__file__).parent / "data" / "h.all.v2023.2.Hs.json.txt"
    if not hallmarks_path.exists():
        raise ValueError(
            "Download the hallmarks genesets json from: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp"
        )
    with open(hallmarks_path) as fd:
        hallmarks = json.load(fd)
    return hallmarks[geneset]["geneSymbols"]
