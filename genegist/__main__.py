from argparse import ArgumentParser
import pickle

from tqdm import tqdm

from genegist.analyzer import Analyzer, Embedding
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
        "--create-dry-run",
        help="Don't actually run the biological process finder, but save the gene summaries to a file",
    )
    parser.add_argument(
        "-a", "--abstracts", help="Also look up abstracts", action="store_true"
    )
    parser.add_argument(
        "-r",
        "--load-dry-run",
        help="Load the gene summaries from a file instead of running the the LLM on them explictly",
    )
    parser.add_argument(
        "--llm",
        help="Specify the LLM to use",
        choices=["gpt-3.5-turbo-1106", "gpt-4-1106-preview"],
        default="gpt-4-1106-preview",
    )
    parser.add_argument(
        "-m",
        "--article",
        help="Get the summary for a given PMID",
    )
    parser.add_argument(
        "-y",
        "--synthetic-generifs",
        help="Create synthetic generifs and save them to a tab-delimited file",
    )
    parser.add_argument(
        "-i",
        "--build-index",
        help="Build an embedding index for all the generifs",
    )

    args = parser.parse_args()

    if args.process:
        generifs = GeneRIFS()
        genes = []
        analyzer = Analyzer(biological_process=args.process, llm=args.llm)
        if args.geneset:
            genes.extend(generifs.get_genes_by_geneset(args.geneset))
        elif args.gene:
            genes.append(args.gene)
        elif args.geneset_file:
            with open(args.geneset_file) as f:
                genes.extend(f.readlines())
        if args.load_dry_run:
            with open(args.load_dry_run, "rb") as f:
                genes = pickle.load(f)
        elif args.create_dry_run:
            with open(args.create_dry_run, "wb") as f:
                summary = analyzer.find_biological_process_from_genes(
                    genes, just_summaries=True, add_abstracts=args.abstracts
                )
                pickle.dump(summary, f)
            return
        print(
            analyzer.find_biological_process_from_genes(
                genes, add_abstracts=args.abstracts
            )
        )
        return

    if args.synthetic_generifs:
        analyzer = Analyzer(llm=args.llm)
        gen = analyzer.create_synthetic_generifs_paired_with_ground_truth(args.gene)
        with open(args.synthetic_generifs, "w") as f:
            f.write("gene id\tgene name\tground truth\tsynth genegist\n")
            for gid, gname, ground, synth in tqdm(
                gen, desc="Writing synthetic generifs", unit=" generif"
            ):
                f.write(f"{gid}\t{gname}\t{ground}\t{synth}\n")

        return

    if args.gene:
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene(args.gene))

    if args.geneset_file:
        generifs = GeneRIFS()
        print(generifs.get_texts_by_gene_set_list(args.geneset_file))

    if args.article:
        analyzer = Analyzer(llm=args.llm)
        print(analyzer.summarize_article(args.article))

    if args.build_index:
        embedding = Embedding()
        embedding.build_index(args.build_index)
