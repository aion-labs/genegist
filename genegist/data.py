from pathlib import Path
from typing import Sequence

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
