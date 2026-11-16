#!/usr/bin/env python3
"""
Mean-pooled ESM-2 embeddings for receptor and precursor sequences.

Frozen model, no fine-tuning: with 386 positive pairs, fine-tuning a 650M
parameter model would memorise the set. The embeddings are features for a small
model trained on top.

Pooling: mean over residue positions, excluding the BOS/EOS tokens. Long
sequences are truncated to --max-len (ESM-2 is trained at 1022 residues); the
report counts how many were affected.

Resumable: writes a .npz per batch into <out>/parts/ and merges at the end.
Re-running skips accessions already embedded.

Usage (Sherlock GPU)
--------------------
    sh_dev                       # or salloc for a longer GPU session
    source $OAK/Tanya/miniconda/bin/activate && conda activate grag
    export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
    nvidia-smi -L                # copy the MIG UUID
    export CUDA_VISIBLE_DEVICES=<MIG UUID>
    python embed_esm.py --batch 4

Usage (CPU)
-----------
    python embed_esm.py --device cpu --batch 2

Outputs:
    data/esm/embeddings.npz      accession -> vector
    data/esm/embed_report.txt
"""

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sequences", default="data/sequences.csv")
    ap.add_argument("--out", default="data/esm")
    ap.add_argument("--model", default="esm2_t33_650M_UR50D",
                    help="esm2_t12_35M_UR50D is a fast smaller alternative")
    ap.add_argument("--batch", type=int, default=4)
    ap.add_argument("--max-len", type=int, default=1022)
    ap.add_argument("--device", default="auto",
                    choices=["auto", "cuda", "cpu"])
    a = ap.parse_args()

    try:
        import torch
        import esm
    except ImportError as e:
        sys.exit(f"missing dependency ({e}). Install with:\n"
                 "  mamba install -c conda-forge pytorch -y\n"
                 "  pip install fair-esm")

    log = []

    def say(m=""):
        print(m, flush=True)
        log.append(m)

    device = a.device
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    say(f"device: {device}")
    if device == "cuda":
        say(f"  gpu: {torch.cuda.get_device_name(0)}")
    say(f"model: {a.model}")

    df = pd.read_csv(a.sequences, dtype=str)
    df["length"] = df["sequence"].str.len()
    say(f"sequences: {len(df)}")

    parts_dir = os.path.join(a.out, "parts")
    os.makedirs(parts_dir, exist_ok=True)

    done = {}
    for fn in sorted(os.listdir(parts_dir)):
        if fn.endswith(".npz"):
            with np.load(os.path.join(parts_dir, fn)) as z:
                for k in z.files:
                    done[k] = z[k]
    if done:
        say(f"resuming: {len(done)} already embedded")

    todo = df[~df["accession"].isin(done)].copy()
    todo = todo.sort_values("length")          # group similar lengths together
    say(f"to embed: {len(todo)}")

    if len(todo):
        model, alphabet = getattr(esm.pretrained, a.model)()
        batch_converter = alphabet.get_batch_converter()
        model = model.to(device).eval()
        n_layers = model.num_layers
        say(f"  layers: {n_layers}, taking the final representation")

        truncated = 0
        t0 = time.time()
        for i in range(0, len(todo), a.batch):
            chunk = todo.iloc[i:i + a.batch]
            data = []
            for _, r in chunk.iterrows():
                seq = r["sequence"][:a.max_len]
                if len(r["sequence"]) > a.max_len:
                    truncated += 1
                data.append((r["accession"], seq))

            _, _, toks = batch_converter(data)
            toks = toks.to(device)
            with torch.no_grad():
                out = model(toks, repr_layers=[n_layers])
            reps = out["representations"][n_layers]

            batch = {}
            for j, (acc, seq) in enumerate(data):
                # exclude BOS (position 0) and EOS (position len+1)
                vec = reps[j, 1:len(seq) + 1].mean(0).float().cpu().numpy()
                batch[acc] = vec
                done[acc] = vec
            np.savez_compressed(
                os.path.join(parts_dir, f"part_{i:05d}.npz"), **batch)

            elapsed = time.time() - t0
            rate = (i + len(chunk)) / max(elapsed, 1e-6)
            left = (len(todo) - i - len(chunk)) / max(rate, 1e-6)
            say(f"  [{i + len(chunk):>4}/{len(todo)}] "
                f"{elapsed / 60:.1f} min elapsed, ~{left / 60:.1f} min left")

        if truncated:
            say(f"  {truncated} sequences truncated to {a.max_len} residues")

    if not done:
        sys.exit("nothing embedded")

    dim = len(next(iter(done.values())))
    say()
    say(f"embeddings: {len(done)} accessions x {dim} dimensions")

    missing = sorted(set(df["accession"]) - set(done))
    if missing:
        say(f"  missing {len(missing)}: {missing[:10]}")

    path = os.path.join(a.out, "embeddings.npz")
    np.savez_compressed(path, **done)
    size = os.path.getsize(path) / 1e6
    say(f"wrote {path} ({size:.1f} MB)")

    # a quick sanity check: nearest neighbours of a known receptor
    accs = sorted(done)
    M = np.stack([done[k] for k in accs])
    M = M / np.linalg.norm(M, axis=1, keepdims=True)
    probe = "P30411" if "P30411" in done else accs[0]
    sims = M @ M[accs.index(probe)]
    order = np.argsort(-sims)[1:6]
    say()
    say(f"sanity check — nearest neighbours of {probe} by cosine similarity:")
    for k in order:
        say(f"  {accs[k]}  {sims[k]:.3f}")
    say("  (for P30411/BDKRB2 the neighbours should be other peptide GPCRs, "
        "not precursors; if precursors rank top, check the pooling)")

    with open(os.path.join(a.out, "embed_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nwrote {a.out}/embed_report.txt")


if __name__ == "__main__":
    main()
