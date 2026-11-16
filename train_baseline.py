#!/usr/bin/env python3
"""
First trained model: can sequence embeddings rank cognate peptide-GPCR pairs?

Evaluated on exactly the splits used for the co-evolution benchmark, so the
numbers are directly comparable (co-evolution reached AUROC ~0.56 under the
primary design).

Design decisions and why
------------------------
* **PCA is fitted on the training folds only.** Fitting it on all pairs would
  leak test-fold structure into the features. The same applies to the scaler.
* **Few features.** 386 positives cannot support 2,560 raw dimensions. Each
  side is reduced to --components, plus their element-wise product as an
  interaction term: 3 * components features.
* **Cold-receptor folds.** Receptors held out whole, as for an orphan.
* **Regularised linear model first.** Gradient boosting is reported alongside,
  but a linear baseline that cannot overfit is the honest reference point.
* **Permutation null**, because a theoretical 0.5 was not trustworthy at this
  sample size earlier in the project.

Usage:
    python train_baseline.py
    python train_baseline.py --components 16 --splits data/splits_class3

Outputs:
    data/model/model_results.csv
    data/model/model_report.txt
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd


def auroc(y, s):
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    keep = ~np.isnan(s)
    y, s = y[keep], s[keep]
    npos, nneg = int((y == 1).sum()), int((y == 0).sum())
    if npos == 0 or nneg == 0:
        return np.nan
    order = np.argsort(s)
    ranks = np.empty(len(s), dtype=float)
    ranks[order] = np.arange(1, len(s) + 1)
    ranks = pd.DataFrame({"s": s, "r": ranks}).groupby("s")["r"] \
             .transform("mean").to_numpy()
    return (ranks[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)


def auprc(y, s):
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    keep = ~np.isnan(s)
    y, s = y[keep], s[keep]
    if (y == 1).sum() == 0:
        return np.nan
    order = np.argsort(-s)
    y = y[order]
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    return float((prec * y).sum() / y.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", default="data/splits_class3")
    ap.add_argument("--pairs", default="data/pairs_resolved.csv")
    ap.add_argument("--embeddings", default="data/esm/embeddings.npz")
    ap.add_argument("--out", default="data/model")
    ap.add_argument("--components", type=int, default=24)
    ap.add_argument("--n-perm", type=int, default=1000)
    a = ap.parse_args()

    try:
        from sklearn.decomposition import PCA
        from sklearn.linear_model import LogisticRegression
        from sklearn.ensemble import HistGradientBoostingClassifier
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import Pipeline
    except ImportError:
        sys.exit("needs scikit-learn:  pip install scikit-learn")

    log = []

    def say(m=""):
        print(m, flush=True)
        log.append(m)

    z = np.load(a.embeddings)
    emb = {k: z[k] for k in z.files}
    dim = len(next(iter(emb.values())))
    say(f"embeddings: {len(emb)} proteins x {dim} dims")

    pairs = pd.read_csv(a.pairs, dtype=str)
    rec_acc = dict(zip(pairs["receptor_gtopdb_id"], pairs["receptor_uniprot"]))
    pep_acc = dict(zip(pairs["peptide_gtopdb_id"], pairs["precursor_uniprot"]))

    files = sorted(glob.glob(os.path.join(a.splits, "splits_*seed*.csv")))
    if not files:
        sys.exit(f"no split files under {a.splits}")
    say(f"split files: {len(files)}  ({a.components} components per side, "
        f"{3 * a.components} features)")

    rows = []
    for path in files:
        m = re.search(r"seed(\d+)", os.path.basename(path))
        seed = int(m.group(1)) if m else 0
        s = pd.read_csv(path, dtype=str)

        # assemble raw feature matrix
        R, P, y, fold, dropped = [], [], [], [], 0
        for _, r in s.iterrows():
            ra = rec_acc.get(r["receptor_gtopdb_id"])
            pa = pep_acc.get(r["peptide_gtopdb_id"])
            if ra not in emb or pa not in emb:
                dropped += 1
                continue
            R.append(emb[ra])
            P.append(emb[pa])
            y.append(int(r["label"]))
            fold.append(int(r["fold"]))
        R, P = np.array(R), np.array(P)
        y, fold = np.array(y), np.array(fold)
        if len(y) == 0:
            say(f"  {os.path.basename(path)}: no usable rows")
            continue

        MODES = ["full", "no_interaction", "receptor_only", "peptide_only"]
        oof = {f"{m}__{mode}": np.full(len(y), np.nan)
               for mode in MODES for m in ("logreg", "gbm")}

        for k in sorted(set(fold)):
            if k < 0:
                continue
            te = fold == k
            tr = (~te) & (fold >= 0)
            if tr.sum() == 0 or len(set(y[tr])) < 2 or te.sum() == 0:
                continue

            # PCA fitted on training rows only -- no leakage
            pr = PCA(n_components=min(a.components, tr.sum() - 1),
                     random_state=seed).fit(R[tr])
            pp = PCA(n_components=min(a.components, tr.sum() - 1),
                     random_state=seed).fit(P[tr])

            def feats(idx, mode="full"):
                r_ = pr.transform(R[idx])
                p_ = pp.transform(P[idx])
                if mode == "receptor_only":
                    # cannot possibly know WHICH peptide; whatever this
                    # scores is the artifact floor
                    return r_
                if mode == "peptide_only":
                    return p_
                if mode == "no_interaction":
                    return np.hstack([r_, p_])
                return np.hstack([r_, p_, r_ * p_])

            for mode in MODES:
                Xtr, Xte = feats(tr, mode), feats(te, mode)

                lr = Pipeline([("sc", StandardScaler()),
                               ("m", LogisticRegression(C=0.1, max_iter=2000,
                                                        random_state=seed))])
                lr.fit(Xtr, y[tr])
                oof[f"logreg__{mode}"][te] = lr.predict_proba(Xte)[:, 1]

                gb = HistGradientBoostingClassifier(
                    max_depth=3, max_iter=200, learning_rate=0.05,
                    l2_regularization=1.0, random_state=seed)
                gb.fit(Xtr, y[tr])
                oof[f"gbm__{mode}"][te] = gb.predict_proba(Xte)[:, 1]

        say()
        say(f"{os.path.basename(path)}  n={len(y)}, "
            f"base rate {y.mean():.3f}, dropped {dropped}")

        rng = np.random.default_rng(seed)
        for name, pred in oof.items():
            ok = ~np.isnan(pred)
            if ok.sum() == 0:
                continue
            yy, ss = y[ok], pred[ok]
            obs = auroc(yy, ss)

            null = np.empty(a.n_perm)
            ff = fold[ok]
            for i in range(a.n_perm):
                yp = yy.copy()
                for k in np.unique(ff):
                    idx = np.where(ff == k)[0]
                    yp[idx] = rng.permutation(yy[idx])
                null[i] = auroc(yp, ss)
            null = null[~np.isnan(null)]
            centre = null.mean()
            lo, hi = np.percentile(null, [2.5, 97.5])
            p = max((np.abs(null - centre) >= abs(obs - centre)).mean(),
                    1.0 / len(null))

            flag = "  <-- exceeds null" if (obs > hi or obs < lo) else ""
            say(f"  {name:<8} AUROC {obs:.3f}  AUPRC {auprc(yy, ss):.3f}  "
                f"null {centre:.3f} [{lo:.3f}–{hi:.3f}]  p={p:.4f}{flag}")
            rows.append({"seed": seed, "model": name, "n": int(ok.sum()),
                         "auroc": obs, "auprc": auprc(yy, ss),
                         "null_centre": centre, "null_lo": lo, "null_hi": hi,
                         "p_value": p})

    if not rows:
        sys.exit("nothing evaluated")

    res = pd.DataFrame(rows)
    os.makedirs(a.out, exist_ok=True)
    res.to_csv(os.path.join(a.out, "model_results.csv"), index=False)

    say()
    say("Mean across seeds")
    for model, g in res.groupby("model"):
        say(f"  {model:<8} AUROC {g.auroc.mean():.3f} "
            f"(seeds {g.auroc.min():.3f}-{g.auroc.max():.3f})  "
            f"AUPRC {g.auprc.mean():.3f}  p {g.p_value.max():.4f} (worst)")

    say()
    say("How to read the ablations. receptor_only and peptide_only are given "
        "ONE side of the pair and cannot know which partner is being scored, "
        "so whatever they achieve is the artifact floor set by the sampling "
        "rather than learned pairing. The full model only demonstrates "
        "pairing to the extent it beats BOTH of them. no_interaction drops "
        "the element-wise product, isolating how much the interaction term "
        "contributes. Co-evolution reached AUROC ~0.56 on these same splits.")

    with open(os.path.join(a.out, "model_report.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print(f"\nwrote {a.out}/model_results.csv and {a.out}/model_report.txt")


if __name__ == "__main__":
    main()
