import json
import os
from pathlib import Path
from typing import Sequence

from Bio import Entrez
import pandas as pd


class GeneRIFS:
    def __init__(self):
        self.df = self.get_generifs()

    def get_generifs(self) -> pd.DataFrame:
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

    def get_texts_by_gene_id(self, gene_id: int) -> Sequence[str]:
        """Get all the texts for the given gene IDs."""
        mask = (self.df["#Tax ID"] == 9606) & (self.df["Gene ID"] == gene_id)
        return self.df[mask]["GeneRIF text"].tolist()

    def get_texts_by_gene_set(self, gene_set: str) -> dict:
        """Get all the texts for the given gene set."""
        genes = []
        gene_ids = []
        for gene in geneset2symbols(gene_set):
            if gene2id(gene) is None:
                raise ValueError(f"Gene {gene} not found.")
            else:
                genes.append(gene)
                gene_ids.append(gene2id(gene))
        return {
            gene: self.get_texts_by_gene_id(gene_id)
            for gene, gene_id in zip(genes, gene_ids)
        }


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
