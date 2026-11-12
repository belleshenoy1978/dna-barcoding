"""Maximum-likelihood phylogeny with IQ-TREE 2 (via the piqtree bindings).

ModelFinder picks the substitution model by BIC unless one is given, and node
support comes from the ultrafast bootstrap (UFBoot). UFBoot values are not read
like standard bootstrap percentages: >= 95 % is the usual "strongly supported"
line, whereas 70 % is the equivalent line for standard bootstrap.

Inference runs in a separate process (see `_iqtree_worker`); on Windows this
avoids the OpenMP runtime clash between IQ-TREE and NumPy/MKL.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

from .config import TREE_DEFAULTS

_WORKER = "barcode_pipeline._iqtree_worker"


def infer_tree(fasta_path: str | Path, model: str | None = None,
               bootstrap_replicates: int | None = None, threads: int = 1,
               seed: int = 1, model_set=None, timeout: int | None = None) -> dict:
    """ModelFinder (optional) + ML tree + UFBoot. Returns a result dictionary.

    Ultrafast bootstrap needs at least 1000 replicates; a smaller value is
    treated as "no bootstrap" by IQ-TREE and reported as such here.
    """
    reps = (TREE_DEFAULTS["bootstrap_replicates"] if bootstrap_replicates is None
            else bootstrap_replicates)
    payload = {"fasta": str(Path(fasta_path).resolve()), "model": model,
               "bootstrap": int(reps or 0), "threads": int(threads), "seed": int(seed),
               "model_set": list(model_set) if model_set else None}
    pkg_parent = str(Path(__file__).resolve().parent.parent)
    env = dict(os.environ, KMP_DUPLICATE_LIB_OK="TRUE")
    env.pop("PYTHONSAFEPATH", None)
    env["PYTHONPATH"] = os.pathsep.join(
        [pkg_parent, *(q for q in env.get("PYTHONPATH", "").split(os.pathsep) if q)])
    proc = subprocess.run([sys.executable, "-m", _WORKER, json.dumps(payload)],
                          capture_output=True, text=True, env=env, timeout=timeout,
                          cwd=pkg_parent)
    if "<<<JSON>>>" not in proc.stdout:
        raise RuntimeError(f"IQ-TREE worker failed (rc={proc.returncode}):\n"
                           f"{proc.stdout[-1500:]}\n{proc.stderr[-1500:]}")
    return json.loads(proc.stdout.split("<<<JSON>>>", 1)[1])


def write_newick(newick: str, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(newick if newick.rstrip().endswith(";") else newick + ";",
                    encoding="utf-8")
    return path


def _read(newick: str):
    import dendropy
    return dendropy.Tree.get(data=newick if newick.rstrip().endswith(";") else newick + ";",
                             schema="newick", preserve_underscores=True)


def midpoint_root_newick(newick: str) -> str:
    """Midpoint-root a Newick string (used when no outgroup is given)."""
    t = _read(newick)
    t.reroot_at_midpoint(update_bipartitions=False)
    return t.as_string(schema="newick", suppress_rooting=True).strip()


def root_with_outgroup(newick: str, outgroup_labels: list[str]) -> str:
    """Root on the MRCA of the given tip labels; falls back to midpoint rooting."""
    t = _read(newick)
    wanted = set(outgroup_labels)
    taxa = [tx for tx in t.taxon_namespace if tx.label in wanted]
    if not taxa:
        return midpoint_root_newick(newick)
    node = t.mrca(taxa=taxa) if len(taxa) > 1 else t.find_node_with_taxon_label(taxa[0].label)
    if node is None or node is t.seed_node or node.edge is None:
        return midpoint_root_newick(newick)
    t.reroot_at_edge(node.edge, update_bipartitions=False)
    return t.as_string(schema="newick", suppress_rooting=True).strip()
