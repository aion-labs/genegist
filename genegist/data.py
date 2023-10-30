from pathlib import Path
import pandas as pd


def get_generifs() -> pd.DataFrame:
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
