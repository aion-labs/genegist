import json
import os
from pathlib import Path
from typing import Sequence

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


def gene2id(gene_name: str) -> int:
    if os.environ.get("NCBI_EMAIL"):
        Entrez.email = os.environ["NCBI_EMAIL"]
    else:
        raise ValueError("Please set the NCBI_EMAIL environment variable.")
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


def get_gene_abstracts(gene_name: str, max_results: int = 100) -> str:
    if os.environ.get("NCBI_EMAIL"):
        Entrez.email = os.environ["NCBI_EMAIL"]
    else:
        raise ValueError("Please set the NCBI_EMAIL environment variable.")
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


def geneset2symbols(geneset: str) -> Sequence[str]:
    hallmarks_path = Path(__file__).parent / "data" / "h.all.v2023.2.Hs.json.txt"
    if not hallmarks_path.exists():
        raise ValueError(
            "Download the hallmarks genesets json from: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp"
        )
    with open(hallmarks_path) as fd:
        hallmarks = json.load(fd)
    return hallmarks[geneset]["geneSymbols"]
