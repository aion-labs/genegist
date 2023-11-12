from typing import Dict

from paperqa import Docs


async def find_biological_process(data: Dict) -> str:
    prompt = "Identify a subset of genes that represents a biological process, "
    "and the role the genes in the process. You want to be as specific as possible."
    "Ideally we want biological pathways."
    index = Docs()
    res = []
    for gene, rif in data.items():
        for i, r in enumerate(rif):
            if r[-1] not in [".", "?", "!"]:
                rif[i] = r + "."
        rif = " ".join(rif)
        if rif.strip() == "":
            continue
        index.add(f"For the gene {gene}: {rif}")

    response = await index.aquery(prompt)
    return response.answer
