# 👓 GeneGist

Researchers often face challenges in deciphering the complex interactions and functions of systems of genes. GeneGist addresses this problem by providing detailed summaries and insights into gene behaviors, interactions, and their roles in biological processes. 

This complexity arises from the vast array of gene interactions, regulatory mechanisms, and the multifaceted roles genes play in biological processes. GeneGist generates in-depth summaries and insights into gene behaviors and interactions, as well as their roles in biological pathways and systems.

GeneGist first scrapes and analyzes academic articles. GeneGist leverages the most advanced Large Language Models (LLMs) available to analyze this information. Using this distilled knowledge it produces biological process summaries. 

GeneGist can also create Gene Reference Into Function (GeneRIFs) directly from scientific literature. GeneRIFs are concise sentence-like annotations, typically written by a human, that describe the function of a gene. GeneGist can construct GeneRIFs using generative AI technology based on LLMs.

## License
Apache License

### Installation

To install GeneGist, ensure you have Python 3.10 or higher. It can be installed via pip:

```bash
pip install genegist
```

## Development 

### Installing Poetry
Poetry is required to handle dependencies and package management. To install Poetry, run:

```bash
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
```

### Setting Up genegist
1. Clone the repository:

   ```bash
   git clone [repository URL]
   cd genegist
   ```

2. Install the dependencies using Poetry:

   ```bash
   poetry install
   ```

## Usage

To use genegist, run the following command:

```bash
poetry run genegist [options]
```

### Options

- `-g GENE`, `--gene GENE`: 
  Look up GeneRIFs for a given gene.

- `-s GENESET`, `--geneset GENESET`: 
  Look up GeneRIFs for a given gene set.

- `-f GENESET_FILE`, `--geneset-file GENESET_FILE`: 
  Look up GeneRIFs for a file containing a list of genes.

- `-p PROCESS`, `--process PROCESS`: 
  Find a biological process for the inputted gene set.

- `-d CREATE_DRY_RUN`, `--create-dry-run CREATE_DRY_RUN`: 
  Don't actually run the biological process finder, but save the gene summaries to a file.

- `-a`, `--abstracts`: 
  Also look up abstracts.

- `-r LOAD_DRY_RUN`, `--load-dry-run LOAD_DRY_RUN`: 
  Load the gene summaries from a file instead of running the LLM on them explicitly.

- `--llm {gpt-3.5-turbo-1106,gpt-4-1106-preview}`: 
  Specify the LLM to use.

- `-m ARTICLE`, `--article ARTICLE`: 
  Get the summary for a given PMID.

- `-t`, `--tasks`: 
  Run a given custom task. Currently only E3 ligase analysis is supported.

- `-y`, `--synthetic-generifs`: 
  Create synthetic generifs and save them to a tab-delimited file.

- `-i`, `--build-index`:
   Build an embedding index for all the generifs.


## Development

### Running Tests
To run tests, use:

```bash
poetry run pytest
```