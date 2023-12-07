from typing import Union, Dict, Iterable

from openai import OpenAI
import tiktoken

from genegist.data import GeneRIFS


def distill(
    full_text: str,
    encoding: tiktoken.Encoding = tiktoken.encoding_for_model("gpt-4-1106-preview"),
    max_tokens: int = 128000,
) -> str:
    """
    Summarizes a given text by breaking it into chunks and processing each chunk using an OpenAI model.

    Args:
        full_text (str): The text to be summarized.
        encoding (tiktoken.Encoding): The encoding to be used for tokenizing the text. Defaults to the encoding suitable for "gpt-4-1106-preview".
        max_tokens (int): The maximum number of tokens that can be processed in a single request. Defaults to 128000.

    Returns:
        str: A summarized version of the provided text.
    """

    prompt = "Summarize the following text:\n\n"
    tokens = encoding.encode(full_text)
    prompt_size = len(encoding.encode(prompt))
    if len(tokens) <= max_tokens:
        return full_text

    chunks = []
    working_chunk = []
    size = 0
    client = OpenAI()

    for i, token in enumerate(tokens):
        working_chunk.append(token)
        size += 1

        if size >= (max_tokens - prompt_size) or i == len(tokens) - 1:
            decoded_chunk = [encoding.decode(part) for part in working_chunk]
            working_chunk = []
            size = 0

            response = client.chat.completions.create(
                messages=[
                    {
                        "role": "user",
                        "content": prompt + decoded_chunk,
                    }
                ],
                model="gpt-3.5-turbo-1106",
            )

            chunks.append(response.choices[0].message.content.strip())

    final_summary = " ".join(chunks)

    return final_summary


def summarize_gene(gene: str, rif: Iterable[str], biological_process: str) -> str:
    """
    Generates a summary for a given gene based on a set of GeneRIFs (Gene Reference Into Function) data.

    Args:
    gene (str): The name of the gene to be summarized.
    rif (Iterable[str]): A collection of GeneRIFs data related to the gene.

    Returns:
    str: A concise, one-paragraph summary of the gene's biological activity and functions.
    """
    client = OpenAI()

    # Ensure each GeneRIF ends with a punctuation mark.
    for i, r in enumerate(rif):
        if r[-1] not in [".", "?", "!"]:
            rif[i] = r + "."

    data = distill(" ".join(rif))

    # Construct the summary prompt.
    summarize_prompt = (
        f"Provide a concise, one-paragraph summary of the biological activity and functions "
        f"of the '{gene}' gene. Base your summary on the following information: "
        f"{data}. Focus on key aspects such as gene expression, regulatory mechanisms, "
        "and its role in cellular processes or disease states, as relevant."
    )

    # Request the summary from the OpenAI API.
    summarize = client.chat.completions.create(
        messages=[{"role": "user", "content": summarize_prompt}],
        model="gpt-4-1106-preview",
    )

    return summarize.choices[0].message.content

def find_biological_process_from_summaries(
    genes: Dict[str, str], biological_process: str
) -> str:
    """
    Identifies genes involved in a specific biological process based on their summaries.

    Args:
    genes (Dict[str, str]): A dictionary where keys are gene names and values are their summaries.
    biological_process (str): The name of the biological process to be investigated.

    Returns:
    str: A narrative describing the subset of genes and their roles in the specified biological process.
    """
    client = OpenAI()

    # Prepare subprompts for each gene.
    gene_subprompts = []
    for gene, summary in genes.items():
        subprompt = f"For the gene {gene}, the biological activity is: {summary}"
        gene_subprompts.append(subprompt)
    
    # Construct the prompt to identify genes involved in the biological process.
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

    # Request the analysis from the OpenAI API.
    bioprocess = client.chat.completions.create(
        messages=[{"role": "user", "content": bioprocess_prompt}],
        model="gpt-4-1106-preview",
    )
    return bioprocess.choices[0].message.content


def find_biological_process_from_genes(
    input_genes: Union[Iterable[Union[str, int]], dict],
    biological_process: str,
    just_summaries: bool = False,
) -> str:
    """
    Determines the involvement of a set of genes in a specific biological process, either by summarizing gene information
    or directly from provided gene summaries.

    Args:
    input_genes (Union[Iterable[Union[str, int]], dict]): A collection of genes or a dictionary of gene summaries.
    biological_process (str): The biological process to investigate.
    just_summaries (bool, optional): If True, returns only the gene summaries. Defaults to False.

    Returns:
    str: Either a dictionary of gene summaries or a narrative about the genes involved in the biological process.
    """
    # Convert the input genes to a dictionary of summaries if not already in dictionary form.
    if isinstance(input_genes, dict):
        genes = input_genes
    else:
        genes = dict()
        for gene in input_genes:
            genes[gene] = summarize_gene(
                gene,
                GeneRIFS().get_texts_by_gene(gene, add_abstracts=True),
            )

    # Return just the summaries if requested.
    if just_summaries:
        return genes

    # Otherwise, find the genes involved in the biological process.
    return find_biological_process_from_summaries(genes, biological_process)
