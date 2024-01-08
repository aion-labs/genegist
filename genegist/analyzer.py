import logging
from time import sleep
import re
from typing import Union, Dict, Iterable, Tuple

from openai import OpenAI, RateLimitError
import numpy as np
from sentence_transformers import SentenceTransformer
import tiktoken
from tqdm import tqdm
import vectorizedb

from genegist.data import (
    GeneRIFS,
    get_gene_abstracts,
    get_article,
    get_abstract,
    id2gene,
    gene2id,
)


class Analyzer:
    def __init__(self, biological_process: str = None, llm: str = "gpt-4-1106-preview"):
        self.biological_process = biological_process
        self.llm = llm

    def call_llm(self, prompt: str, llm: str) -> str:
        """
        Calls an OpenAI language model to process a given prompt.

        Args:
            prompt (str): The prompt to be processed.
            llm (str): The language model to be used.

        Returns:
            str: The processed prompt.
        """

        client = OpenAI()

        try:
            # Request the summary from the OpenAI API.
            result = client.chat.completions.create(
                messages=[{"role": "user", "content": prompt}], model=self.llm
            )

        except RateLimitError as e:
            regex = re.compile(r"Please try again in (\d+)s.")
            match = regex.search(e.body["error"]["message"])
            wait_time = int(match.group(1) + 5)
            logging.warning(f"Rate limit reached, waiting {wait_time} seconds.")
            if match:
                sleep(wait_time)
                result = client.chat.completions.create(
                    messages=[{"role": "user", "content": prompt}],
                    model=llm,
                )

        return result.choices[0].message.content.strip()

    def summarize_article(self, pmid: str) -> str:
        """
        Summarizes a given article.

        Args:
            pmid (str): The PubMed ID of the article to be summarized.

        Returns:
            str: A summary of the article.
        """

        article = get_article(pmid)
        if article is None:
            return "No article found."

        article = self.distill(article)

        prompt = f"Summarize the following article:\n\n{article}"

        return self.call_llm(prompt, self.llm)

    def distill(
        self,
        full_text: str,
        encoding: tiktoken.Encoding = None,
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

        if encoding is None:
            encoding = tiktoken.encoding_for_model(self.llm)
        if self.llm == "gpt-4-1106-preview":
            max_tokens = 128000
        else:
            max_tokens = 16000
        gpt35_max_tokens = 16000

        prompt = "Summarize the following text:\n\n"
        tokens = encoding.encode(full_text)
        prompt_size = len(encoding.encode(prompt))
        if len(tokens) <= max_tokens:
            logging.debug(
                "distill: Text is short enough to be summarized in one request."
            )
            return full_text

        chunks = []
        working_chunk = []
        size = 0

        for i, token in enumerate(tokens):
            working_chunk.append(token)
            size += 1

            if size >= (gpt35_max_tokens - prompt_size) or i == len(tokens) - 1:
                if i == len(tokens) - 1:
                    logging.debug(f"distill: Reached end of text, chunk size {size}.")
                else:
                    logging.debug(f"distill: Working chunk size {size}.")

                decoded_chunk = encoding.decode(working_chunk)
                working_chunk = []
                size = 0

                chunks.append(
                    self.call_llm(prompt + decoded_chunk, "gpt-3.5-turbo-1106")
                )

        final_summary = " ".join(chunks)

        return final_summary

    def create_generifs_from_abstracts(
        self, gene: str, abstracts: str
    ) -> Iterable[str]:
        """
        Creates GeneRIFs (Gene Reference Into Function) data from a set of abstracts.

        Args:
        gene (str): The name of the gene to be summarized.
        abstracts (str): A collection of abstracts related to the gene.

        Returns:
        Iterable[str]: A collection of GeneRIFs data related to the gene.
        """

        abstracts = self.distill(abstracts)

        prompt = (
            f"Provide exactly one concise, very short one-sentence summary of the biological activity and functions "
            f"of the '{gene}' gene related to {self.biological_process}, similar to a GeneRIF. Base your "
            "summary on the following information: "
            f"{abstracts}. Focus on key aspects such as gene expression, regulatory mechanisms, "
            "and its role in cellular processes or disease states, as relevant. Keep it as short as possible, like the size of a tweet."
        )

        result = self.call_llm(prompt, self.llm).split(".")
        result = [r.strip() for r in result if r != ""]
        return result[0]

    def summarize_gene(self, gene: str, rif: Iterable[str]) -> str:
        """
        Generates a summary for a given gene based on a set of GeneRIFs (Gene Reference Into Function) data.

        Args:
        gene (str): The name of the gene to be summarized.
        rif (Iterable[str]): A collection of GeneRIFs data related to the gene.

        Returns:
        str: A concise, one-paragraph summary of the gene's biological activity and functions.
        """

        # Ensure each GeneRIF ends with a punctuation mark.
        for i, r in enumerate(rif):
            if r[-1] not in [".", "?", "!"]:
                rif[i] = r + "."

        data = self.distill(" ".join(rif))

        # Construct the summary prompt.
        summarize_prompt = (
            f"Provide a concise, one-paragraph summary of the biological activity and functions "
            f"of the '{gene}' gene related to f{self.biological_process}. Base your "
            "summary on the following information: "
            f"{data}. Focus on key aspects such as gene expression, regulatory mechanisms, "
            "and its role in cellular processes or disease states, as relevant."
        )

        return self.call_llm(summarize_prompt, self.llm)

    def find_biological_process_from_summaries(self, genes: Dict[str, str]) -> str:
        """
        Identifies genes involved in a specific biological process based on their summaries.

        Args:
        genes (Dict[str, str]): A dictionary where keys are gene names and values are their summaries.
        biological_process (str): The name of the biological process to be investigated.

        Returns:
        str: A narrative describing the subset of genes and their roles in the specified biological process.
        """

        # Prepare subprompts for each gene.
        gene_subprompts = []
        for gene, summary in genes.items():
            subprompt = f"For the gene {gene}, the biological activity is: {summary}"
            gene_subprompts.append(subprompt)

        # TODO: Add a subprompt connecting different biological process.

        # Construct the prompt to identify genes involved in the biological process.
        bioprocess_prompt = (
            "Please identify and describe the subset of genes specifically involved in the "
            f"'{self.biological_process}' biological process, using the provided information. "
            "For each gene, detail its unique role and function within this process. "
            f"Exclude any genes not directly related to '{self.biological_process}'. "
            "Your response should be structured as a clear, comprehensive narrative, "
            "detailing the biological pathway in a prose format. Aim for precision and "
            "specificity in your explanation. Provide a mechanistic hypothesis for the role of each gene in "
            f"{self.biological_process}"
            f"{f'{chr(92)}n'.join(gene_subprompts)}"
        )

        # Request the analysis from the OpenAI API.
        return self.call_llm(bioprocess_prompt, self.llm)

    def create_synthetic_generifs_paired_with_ground_truth(
        self, gene: str = None
    ) -> Iterable[Tuple[str, str, str, str]]:
        """
        Creates synthetic GeneRIFs (Gene Reference Into Function) data from a set of abstracts.

        Args:
            gene (str): The name of the gene to be summarized. Default is all genes.

        Returns:
            Iterable[Tuple[str, str, str, str]]: A collection of tuples of structure: NCBI gene id,
            gene name, ground truth GeneRIF, synthetic GeneRIF.
        """
        generifs = GeneRIFS()
        df = generifs.get_generifs()

        df = df[df["#Tax ID"] == 9606]  # Only human genes
        df = df[df["PubMed ID (PMID) list"].str.len() > 0]  # Only genes with PMIDs
        df = df[df["GeneRIF text"].str.len() > 0]  # Only genes with GeneRIFs
        if gene is not None:
            df = df[df["Gene ID"] == gene2id(gene)]

        for _, row in df.iterrows():
            abstract = get_abstract(row["PubMed ID (PMID) list"])
            gene_id = row["Gene ID"]
            gene_name = id2gene(gene_id)
            if abstract is None or gene_name is None:
                continue
            yield gene_id, gene_name, row[
                "GeneRIF text"
            ], self.create_generifs_from_abstracts(gene_name, abstract)

    def find_biological_process_from_genes(
        self,
        input_genes: Union[Iterable[Union[str, int]], dict],
        just_summaries: bool = False,
        add_abstracts: bool = False,
    ) -> Union[str, dict]:
        """
        Determines the involvement of a set of genes in a specific biological process, either by summarizing gene information
        or directly from provided gene summaries.

        Args:
        input_genes (Union[Iterable[Union[str, int]], dict]): A collection of genes or a dictionary of gene summaries.
        biological_process (str): The biological process to investigate.
        just_summaries (bool, optional): If True, returns only the gene summaries. Defaults to False.
        add_abstracts (bool, optional): If True, also uses abstracts to create synthetic GeneRIFs. Defaults to False.

        Returns:
        Union[str, dict]: Either a dictionary of gene summaries or a narrative about the genes involved in the biological process.
        """
        # Convert the input genes to a dictionary of summaries if not already in dictionary form.
        if isinstance(input_genes, dict):
            genes = input_genes
        else:
            genes = dict()
            for gene in input_genes:
                rifs = GeneRIFS().get_texts_by_gene(gene)
                if add_abstracts:
                    abstracts = get_gene_abstracts(gene)
                    rifs.extend(self.create_generifs_from_abstracts(gene, abstracts))
                genes[gene] = self.summarize_gene(
                    gene,
                    rifs,
                )

        # Return just the summaries if requested.
        if just_summaries:
            return genes

        # Otherwise, find the genes involved in the biological process.
        return self.find_biological_process_from_summaries(genes)


class Embedding:
    def __init__(self):
        self.embedding = SentenceTransformer("all-MiniLM-L6-v2")
        generifs = GeneRIFS().get_generifs()
        self.generifs = generifs[generifs["#Tax ID"] == 9606]

    def get_embedding(self, gene: Union[str, int]) -> np.ndarray:
        """
        Get the embedding for the given gene.

        Args:
            gene (Union[str, int]): The gene to be embedded by gene ID or gene name.

        Returns:
            np.ndarray: The embedding of the gene.
        """
        if isinstance(gene, str):
            gene = gene2id(gene)

        texts = self.generifs[self.generifs["Gene ID"] == gene]["GeneRIF text"].tolist()
        if len(texts) == 0:
            return None
        return self.embedding.encode(texts)

    def build_index(self, index_path: str) -> None:
        """
        Build the vector space for all the GeneRIFs.

        Args:
            index_path (str): The path to save the index.
        """

        db = vectorizedb.Database(index_path, dim=384, max_elements=9999999)

        try:
            for gene in tqdm(
                self.generifs["Gene ID"].unique(), desc="Building index", unit=" gene"
            ):
                texts = self.generifs[self.generifs["Gene ID"] == gene][
                    "GeneRIF text"
                ].tolist()
                gene_name = id2gene(gene)

                for idx, text in enumerate(texts, 1):
                    embedding = self.embedding.encode(text)
                    if embedding is not None:
                        db[f"{gene}-gt-{idx}"] = (
                            embedding,
                            {
                                "gene_id": int(gene),
                                "gene_name": gene_name,
                                "type": "ground_truth",
                                "text": text,
                            },
                        )
        except KeyboardInterrupt:
            print("Stopping index creation early.")
