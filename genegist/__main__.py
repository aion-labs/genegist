from argparse import ArgumentParser

from genegist.data import GeneRIFS, gene2id


def main():
    parser = ArgumentParser()
    parser.add_argument("-g", "--gene", help="Look up GeneRIFS for a given gene")
    parser.add_argument("-s", "--geneset", help="Look up GeneRIFS for a given gene set")
    parser.add_argument(
        "-a", "--abstracts", help="Also look up abstracts", action="store_true"
    )
    args = parser.parse_args()

    if args.gene:
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene(args.gene, args.abstracts))

    if args.geneset:
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene_set(args.geneset, args.abstracts))
