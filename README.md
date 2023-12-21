# 👓 GeneGist

"genegist" is a Python package designed for generating summaries of groups of genes. This tool leverages LLMs to analyze gene-related data effectively.

## License
Apache License

## Installation

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

## Development

### Running Tests
To run tests, use:

```bash
poetry run pytest
```