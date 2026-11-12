"""Worker process that runs IQ-TREE 2 (piqtree) and returns JSON on stdout.

Kept in its own process for two reasons: IQ-TREE's bundled OpenMP runtime
clashes with the one NumPy/MKL loads on Windows (KMP_DUPLICATE_LIB_OK is set
here, before any import), and an IQ-TREE abort then cannot take the calling
session down with it.
"""
import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import json  # noqa: E402
import sys  # noqa: E402


def newick_with_support(tree) -> str:
    """Newick string whose internal node labels are the UFBoot percentages."""
    nwk = tree.get_newick(with_distances=True, with_node_names=True)
    for node in tree.postorder():
        if node.is_tip() or not node.name:
            continue
        sup = getattr(node, "support", None)
        label = f"{sup:g}" if sup is not None else ""
        nwk = nwk.replace(f"){node.name}:", f"){label}:")
        nwk = nwk.replace(f"){node.name};", f"){label};")
    return nwk


def main() -> None:
    args = json.loads(sys.argv[1])
    import cogent3
    import piqtree

    aln = cogent3.load_aligned_seqs(args["fasta"], moltype="dna")
    threads = args.get("threads") or 1
    seed = args.get("seed", 1)

    mf_out = {}
    model = args.get("model")
    if not model:
        ms = args.get("model_set")
        mf = piqtree.model_finder(aln, model_set=set(ms) if ms else None,
                                  num_threads=threads, rand_seed=seed)
        mf_out = {"best_aic": str(mf.best_aic), "best_aicc": str(mf.best_aicc),
                  "best_bic": str(mf.best_bic)}
        model = mf_out["best_bic"]

    reps = args.get("bootstrap") or 0
    reps = reps if reps >= 1000 else 0
    tree = piqtree.build_tree(aln, model, rand_seed=seed,
                              bootstrap_replicates=(reps or None),
                              num_threads=threads)
    params = dict(getattr(tree, "params", {}) or {})
    out = {
        "newick": newick_with_support(tree),
        "model": str(model),
        "model_finder": mf_out,
        "bootstrap_replicates": reps or None,
        "n_tips": len(tree.get_tip_names()),
        "log_likelihood": params.get("lnL"),
        "support_type": "UFBoot (%)" if reps else "none",
    }
    sys.stdout.write("<<<JSON>>>" + json.dumps(out, default=str))


if __name__ == "__main__":
    main()
