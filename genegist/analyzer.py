from typing import Dict
import io
from paperqa import Docs


async def find_biological_process(data: Dict) -> str:
    prompt = (
        "Identify a subset of genes that represents a biological process, "
        "and the role the genes in the process. You want to be as specific as possible. "
        "Ideally, we want biological pathways."
    )
    index = Docs()

    for gene, rif in data.items():
        # Process and clean up the text for each gene
        for i, r in enumerate(rif):
            if r[-1] not in [".", "?", "!"]:
                rif[i] = r + "."
        rif = " ".join(rif)
        if rif.strip() == "":
            continue

        # Create an in-memory bytes-based file-like object for each gene's information
        with io.BytesIO(f"For the gene {gene}: {rif}".encode()) as gene_file:
            index.add_file(gene_file, dockey=gene)

    response = await index.aquery(prompt)
    return response.answer
