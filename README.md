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
- `-g` or `--gene`: Look up GeneRIFS for a given gene
- `-s` or `--geneset`: Look up GeneRIFS for a given gene set
- `-f` or `--geneset-file`: Look up GeneRIFS for file containing a list of genes
- `-p` or `--process`: Find a biological process for the inputted gene set
- `-d` or `--dry-run`: Print the texts without running the biological process finder
- `-a` or `--abstracts`: Also look up abstracts

## Development

### Running Tests
To run tests, use:

```bash
poetry run pytest
```