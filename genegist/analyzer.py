from typing import Union, Dict
import io
from uuid import uuid4
from paperqa import Docs


async def find_biological_process(data: Union[Dict, str], dry_run: bool = False) -> str:
    prompt = (
        "Identify a subset of genes that represents a biological process, "
        "and the role the genes in the process. You want to be as specific as possible. "
        "Ideally, we want biological pathways."
    )
    index = Docs()

    if dry_run:
        index = []

    if isinstance(data, str):
        data = data.splitlines()
        for line in data:
            if line.strip() == "":
                continue
            with io.BytesIO(line.encode()) as gene_file:
                index.add_file(gene_file, dockey=str(uuid4()))
    else:
        for gene, rif in data.items():
            # Process and clean up the text for each gene
            for i, r in enumerate(rif):
                if r[-1] not in [".", "?", "!"]:
                    rif[i] = r + "."
            rif = " ".join(rif)
            if rif.strip() == "":
                continue

            text = f"For the gene {gene}: {rif}"
            if dry_run:
                index.append(text)
                continue
            # Create an in-memory bytes-based file-like object for each gene's information
            with io.BytesIO(text.encode()) as gene_file:
                index.add_file(gene_file, dockey=gene)

    if dry_run:
        return "\n".join(index)

    response = await index.aquery(prompt)
    return response.answer
