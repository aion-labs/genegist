from argparse import ArgumentParser

from genegist.data import GeneRIFS, gene2id


def main():
    parser = ArgumentParser()
    parser.add_argument("-g", "--gene", help="Look up GeneGists for a given gene")
    args = parser.parse_args()

    if args.gene:
        gene_id = gene2id(args.gene)
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene_id(gene_id))
