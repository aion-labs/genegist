from argparse import ArgumentParser
import asyncio

from genegist.analyzer import find_biological_process
from genegist.data import GeneRIFS, gene2id


async def async_main():
    parser = ArgumentParser()
    parser.add_argument("-g", "--gene", help="Look up GeneRIFS for a given gene")
    parser.add_argument("-s", "--geneset", help="Look up GeneRIFS for a given gene set")
    parser.add_argument(
        "-a", "--abstracts", help="Also look up abstracts", action="store_true"
    )
    parser.add_argument(
        "-p", "--process", help="Find a biological process", action="store_true"
    )

    args = parser.parse_args()

    if args.gene:
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene(args.gene, args.abstracts))

    if args.geneset:
        generifs = GeneRIFS()
        texts = generifs.get_texts_by_gene_set(args.geneset, args.abstracts)
        if not args.process:
            print(texts)
            return
        else:
            print(await find_biological_process(texts))


def main():
    asyncio.run(async_main())
