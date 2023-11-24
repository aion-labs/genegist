from argparse import ArgumentParser

from genegist.analyzer import find_biological_process_from_genes
from genegist.data import GeneRIFS


def main():
    parser = ArgumentParser()
    parser.add_argument("-g", "--gene", help="Look up GeneRIFS for a given gene")
    parser.add_argument("-s", "--geneset", help="Look up GeneRIFS for a given gene set")
    parser.add_argument(
        "-f",
        "--geneset-file",
        help="Look up GeneRIFS for file containing a list of genes",
    )
    parser.add_argument(
        "-p",
        "--process",
        help="Find a biological process for the inputed gene set",
    )
    parser.add_argument(
        "-d",
        "--dry-run",
        help="Don't actually run the biological process finder, but print the texts",
        action="store_true",
    )
    parser.add_argument(
        "-a", "--abstracts", help="Also look up abstracts", action="store_true"
    )
    parser.add_argument(
        "-r",
        "--dry-run-file",
        help="Use a dry run file to find biological processes",
    )

    args = parser.parse_args()

    if args.dry_run_file or args.dry_run or args.geneset:
        raise RuntimeError("Some features are temporarily unavailable.")

    if args.process:
        generifs = GeneRIFS()
        genes = []
        if args.gene:
            genes.append(args.gene)
        if args.geneset_file:
            with open(args.geneset_file) as f:
                genes.extend(f.readlines())
        print(find_biological_process_from_genes(genes, args.process))
        return

    if args.gene:
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene(args.gene, args.abstracts))

    if args.geneset_file:
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene_set_list(args.geneset_file, args.abstracts))
