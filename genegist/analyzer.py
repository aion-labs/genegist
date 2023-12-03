from typing import Union, Dict, Iterable

from openai import OpenAI

from genegist.data import GeneRIFS


def summarize_gene(gene: str, rif: Iterable[str]) -> str:
    client = OpenAI()

    for i, r in enumerate(rif):
        if r[-1] not in [".", "?", "!"]:
            rif[i] = r + "."

    # Hacky way to limit the number of tokens
    # TODO: Switch to GPT-4 when OpenAI allows us to use it
    rif = rif[:430]

    summerize_prompt = f"In concise paragraph, summarize the biological activity of the {gene} gene from the following GeneRIFS: {' '.join(rif)}"
    summerize = client.chat.completions.create(
        messages=[{"role": "user", "content": summerize_prompt}],
        model="gpt-4-1106-preview",
    )

    return summerize.choices[0].message.content


def find_biological_process_from_summaries(
    genes: Dict[str, str], biological_process: str
) -> str:
    client = OpenAI()

    gene_subprompts = []
    for gene, summary in genes.items():
        subprompt = f"For the gene {gene}, the biological activity is: {summary}"
        gene_subprompts.append(subprompt)
    
    bioprocess_prompt = (
        "Please identify and describe the subset of genes specifically involved in the "
        f"'{biological_process}' biological process, using the provided information. "
        "For each gene, detail its unique role and function within this process. "
        f"Exclude any genes not directly related to '{biological_process}'. "
        "Your response should be structured as a clear, comprehensive narrative, "
        "detailing the biological pathway in a prose format. Aim for precision and "
        "thoroughness in your explanation."
        f"{f'{chr(92)}n'.join(gene_subprompts)}"
    )

    bioprocess = client.chat.completions.create(
        messages=[{"role": "user", "content": bioprocess_prompt}],
        model="gpt-4-1106-preview",
    )
    return bioprocess.choices[0].message.content


def find_biological_process_from_genes(
    input_genes: Union[Iterable[Union[str, int]], dict], biological_process: str, just_summaries: bool = False
) -> str:
    if isinstance(input_genes, dict):
        genes = input_genes
    else:
        genes = dict()
        for gene in input_genes:
            genes[gene] = summarize_gene(gene, GeneRIFS().get_texts_by_gene(gene))
    if just_summaries:
        return genes
    return find_biological_process_from_summaries(genes, biological_process)
