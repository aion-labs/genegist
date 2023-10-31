from argparse import ArgumentParser

from genegist.data import GeneRIFS, gene2id


def main():
    parser = ArgumentParser()
    parser.add_argument("-g", "--gene", help="Look up GeneRIFS for a given gene")
    parser.add_argument("-s", "--geneset", help="Look up GeneRIFS for a given gene set")
    args = parser.parse_args()

    if args.gene:
        gene_id = gene2id(args.gene)
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene_id(gene_id))

    if args.geneset:
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene_set(args.geneset))
