"""
ESM-2 embedding extractor -- run ONCE, cache to disk.

This is the only script that needs torch + fair-esm. It runs each peptide
through a FROZEN ESM-2 model (weights never change -- this is feature
extraction / transfer learning, NOT fine-tuning) and mean-pools the per-residue
embeddings into one fixed-length vector per peptide. Results are cached to
results/esm_embeddings.npz so all downstream modeling reads the cache and never
needs torch again.

Model choice: esm2_t12_35M_UR50D (35M params). Small and fast on CPU, and
appropriate for a <1000-peptide dataset -- a bigger model's richer embeddings
tend to OVERFIT at this data scale, so small is the honest choice, not a
compromise.

Install (once):
    pip install torch fair-esm

Run (once, ~few minutes on CPU):
    python3 src/embed_esm.py pepic50.txt
Also embed your own peptides so they can be predicted later:
    python3 src/embed_esm.py --seqs WYNAW NAW

Saves: results/esm_embeddings.npz  (dict: sequence -> vector)
"""

from __future__ import annotations
import argparse
import os
import numpy as np


def load_cache(path):
    if os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        return {k: v for k, v in zip(d["seqs"], d["vecs"])}
    return {}


def save_cache(path, cache):
    seqs = list(cache.keys())
    vecs = np.stack([cache[s] for s in seqs]) if seqs else np.zeros((0, 0))
    np.savez_compressed(path, seqs=np.array(seqs, dtype=object), vecs=vecs)


def embed_sequences(seqs, model_name="esm2_t12_35M_UR50D", batch_size=16):
    """Mean-pooled ESM-2 embeddings. Imports torch/esm lazily so the rest of
    the project never depends on them."""
    import torch
    import esm

    print(f"loading {model_name} (first run downloads ~150MB)...")
    model, alphabet = getattr(esm.pretrained, model_name)()
    bc = alphabet.get_batch_converter()
    model.eval()
    # figure out which layer index to read (last layer)
    n_layers = model.num_layers

    out = {}
    for i in range(0, len(seqs), batch_size):
        chunk = seqs[i:i + batch_size]
        data = [(str(j), s) for j, s in enumerate(chunk)]
        _, _, toks = bc(data)
        with torch.no_grad():
            rep = model(toks, repr_layers=[n_layers])["representations"][n_layers]
        for k, s in enumerate(chunk):
            # mean-pool over real residues (skip BOS at 0 and EOS at len+1)
            vec = rep[k, 1:len(s) + 1].mean(0).cpu().numpy()
            out[s] = vec.astype(np.float32)
        print(f"  embedded {min(i+batch_size, len(seqs))}/{len(seqs)}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datafile", nargs="?", default=None,
                    help="AHTPDB file to embed all peptides from")
    ap.add_argument("--seqs", nargs="+", default=None,
                    help="explicit sequences to embed (e.g. your peptides)")
    ap.add_argument("--cache", default="results/esm_embeddings.npz")
    ap.add_argument("--model", default="esm2_t12_35M_UR50D")
    args = ap.parse_args()

    os.makedirs("results", exist_ok=True)
    cache = load_cache(args.cache)
    print(f"cache has {len(cache)} embeddings already")

    # collect sequences to embed
    need = []
    if args.datafile:
        from load_ahtpdb import read_ahtpdb_file, clean_ahtpdb
        df = read_ahtpdb_file(args.datafile)
        clean, _ = clean_ahtpdb(df, seq_col="seq", ic50_col="ic50",
                                sbp_col="bp", verbose=False)
        need += list(clean["sequence"])
    if args.seqs:
        need += [s.strip().upper() for s in args.seqs]

    todo = sorted(set(s for s in need if s not in cache))
    print(f"{len(todo)} new sequences to embed")
    if todo:
        new = embed_sequences(todo, model_name=args.model)
        cache.update(new)
        save_cache(args.cache, cache)
        print(f"saved {len(cache)} total embeddings to {args.cache}")
    else:
        print("nothing new to embed; cache is up to date")


if __name__ == "__main__":
    main()
