from argparse import ArgumentParser
import asyncio

from genegist.analyzer import find_biological_process
from genegist.data import GeneRIFS, gene2id


async def async_main():
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
        action="store_true",
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

    args = parser.parse_args()

    if args.gene:
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene(args.gene, args.abstracts))

    if args.geneset or args.geneset_file:
        generifs = GeneRIFS()
        if args.geneset_file:
            with open(args.geneset_file) as f:
                genes = f.read().splitlines()
            texts = generifs.get_texts_by_gene_set_list(genes, args.abstracts)
        else:
            texts = generifs.get_texts_by_gene_set(args.geneset, args.abstracts)
        if not args.process:
            print(texts)
            return
        else:
            print(await find_biological_process(texts, args.dry_run))


def main():
    asyncio.run(async_main())
