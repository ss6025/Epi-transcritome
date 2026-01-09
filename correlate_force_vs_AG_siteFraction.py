#!/usr/bin/env python3
import os
import argparse
import math
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from intervaltree import IntervalTree
from collections import defaultdict


def build_force_intervaltrees(force_tsv: str):
    """
    Forces TSV required columns:
      chr, start, end, max_mean_dsRNA_force

    Assumption:
      - start is 0-based, end is 1-based exclusive (BED-like).
    """
    df = pd.read_csv(force_tsv, sep="\t", dtype={"chr": str})

    required = {"chr", "start", "end", "max_mean_dsRNA_force"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Forces TSV missing columns: {sorted(missing)}")

    df["start"] = pd.to_numeric(df["start"], errors="raise")
    df["end"] = pd.to_numeric(df["end"], errors="raise")
    df["max_mean_dsRNA_force"] = pd.to_numeric(df["max_mean_dsRNA_force"], errors="coerce")

    trees = {}
    meta = {}

    def make_key(row):
        repid = row.get("rep.id", "NA")
        strand = row.get("strand", ".")
        return f"{row['chr']}:{int(row['start'])}-{int(row['end'])}|{repid}|{strand}"

    for _, row in df.iterrows():
        chrom = row["chr"]
        s = int(row["start"])
        e = int(row["end"])
        key = make_key(row)

        if chrom not in trees:
            trees[chrom] = IntervalTree()
        trees[chrom].addi(s, e, key)
        meta[key] = row.to_dict()

    return trees, meta


def allsubs_has_AG(allsubs: str) -> bool:
    """
    True if AllSubs contains an AG token.
    Accepts formats like: "AG" or "AG,TC" or "AG;TC" etc.
    """
    if not isinstance(allsubs, str):
        return False
    s = allsubs.strip()
    if s == "" or s == "-":
        return False
    toks = []
    for chunk in s.replace(";", ",").replace("|", ",").split(","):
        c = chunk.strip()
        if c:
            toks.append(c)
    return ("AG" in toks) or (s == "AG")


def is_any_edit(allsubs: str) -> bool:
    """Any edited site = AllSubs not '-' and not empty."""
    if not isinstance(allsubs, str):
        return False
    s = allsubs.strip()
    return (s != "") and (s != "-")


def correlate(x, y):
    """Return (pearson_r, spearman_r, n) using pandas corr."""
    s = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(s) < 3:
        return (math.nan, math.nan, len(s))
    pear = s["x"].corr(s["y"], method="pearson")
    spear = s["x"].corr(s["y"], method="spearman")
    return (pear, spear, len(s))


def safe_qcut(series, q, labels=None):
    """qcut with fallback to cut for tied values."""
    try:
        return pd.qcut(series, q=q, labels=labels, duplicates="drop")
    except Exception:
        return pd.cut(series, bins=q, labels=labels)


def main():
    ap = argparse.ArgumentParser(
        description="Correlate max_mean_dsRNA_force with fraction of edited sites that are AG (AllSubs-based, not read-based)."
    )
    ap.add_argument("--forces", required=True, help="Forces TSV (repeat intervals + max_mean_dsRNA_force)")
    ap.add_argument("--reditools_bed", required=True, help="REDItools .bed (tab-delimited) with header")
    ap.add_argument("--out_prefix", required=True, help="Output prefix (no extension)")

    ap.add_argument("--min_cov", type=int, default=1,
                    help="Minimum Coverage-q30 to keep a site (default 1). Set 0 to disable.")
    ap.add_argument("--min_edited_sites_per_interval", type=int, default=10,
                    help="Minimum # edited sites in interval to keep it (default 10)")
    ap.add_argument("--bins", type=int, default=0,
                    help="Number of force bins for curve (0=auto: 30 if >=300 else 15)")
    ap.add_argument("--max_rows", type=int, default=0,
                    help="Debug: stop after N rows after filtering (0=no limit)")

    args = ap.parse_args()

    trees, meta = build_force_intervaltrees(args.forces)

    # Aggregate per interval:
    #   edited_sites = AllSubs != '-'
    #   ag_sites = AllSubs contains AG
    agg = defaultdict(lambda: {
        "edited_sites": 0,
        "ag_sites": 0,
    })

    bed = pd.read_csv(args.reditools_bed, sep="\t", low_memory=True)

    needed = {"Region", "Position", "AllSubs"}
    missing = needed - set(bed.columns)
    if missing:
        raise ValueError(f"REDItools BED missing columns: {sorted(missing)}")

    # Optional coverage filter if present
    if "Coverage-q30" in bed.columns:
        bed["Coverage-q30"] = pd.to_numeric(bed["Coverage-q30"], errors="coerce")
        if args.min_cov > 0:
            bed = bed[bed["Coverage-q30"].fillna(0) >= args.min_cov]

    bed["Position"] = pd.to_numeric(bed["Position"], errors="coerce")
    bed = bed[bed["Position"].notna()]

    # Keep only edited sites (denominator)
    bed["AllSubs"] = bed["AllSubs"].astype(str)
    bed = bed[bed["AllSubs"].apply(is_any_edit)]

    if args.max_rows and args.max_rows > 0:
        bed = bed.head(args.max_rows)

    for _, r in bed.iterrows():
        chrom = str(r["Region"])
        pos1 = int(r["Position"])
        pos0 = pos1 - 1  # assume REDItools Position is 1-based

        if chrom not in trees:
            continue

        hits = trees[chrom].at(pos0)
        if not hits:
            continue

        allsubs = r["AllSubs"]
        is_ag = allsubs_has_AG(allsubs)

        for iv in hits:
            key = iv.data
            agg[key]["edited_sites"] += 1
            if is_ag:
                agg[key]["ag_sites"] += 1

    # Build per-interval table
    rows = []
    for key, a in agg.items():
        if a["edited_sites"] == 0:
            continue
        m = meta.get(key, {})
        force = m.get("max_mean_dsRNA_force", math.nan)
        ag_frac = a["ag_sites"] / a["edited_sites"] if a["edited_sites"] > 0 else math.nan

        rows.append({
            "interval_key": key,
            "chr": m.get("chr", key.split(":")[0]),
            "start": int(m.get("start", -1)) if not pd.isna(m.get("start", np.nan)) else -1,
            "end": int(m.get("end", -1)) if not pd.isna(m.get("end", np.nan)) else -1,
            "rep.id": m.get("rep.id", "NA"),
            "Family": m.get("Family", m.get("family", "NA")),
            "Class": m.get("Class", m.get("class", "NA")),
            "strand": m.get("strand", "."),
            "length": m.get("length", np.nan),
            "max_mean_dsRNA_force": force,
            "edited_sites_in_interval": a["edited_sites"],
            "AG_sites_in_interval": a["ag_sites"],
            "AG_edit_fraction": ag_frac,  # <-- what you want
        })

    outdf = pd.DataFrame(rows)

    # Filter for correlation / plots
    keep = outdf.copy()
    keep["max_mean_dsRNA_force"] = pd.to_numeric(keep["max_mean_dsRNA_force"], errors="coerce")
    keep["AG_edit_fraction"] = pd.to_numeric(keep["AG_edit_fraction"], errors="coerce")

    keep = keep[
        (keep["edited_sites_in_interval"] >= args.min_edited_sites_per_interval) &
        keep["max_mean_dsRNA_force"].notna() &
        keep["AG_edit_fraction"].notna()
    ].copy()

    pear, spear, n = correlate(keep["max_mean_dsRNA_force"], keep["AG_edit_fraction"])

    # Write outputs
    out_dir = os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    out_tsv = args.out_prefix + ".per_interval.tsv"
    keep_tsv = args.out_prefix + ".per_interval.filtered.tsv"
    stats_tsv = args.out_prefix + ".stats.tsv"

    outdf.to_csv(out_tsv, sep="\t", index=False)
    keep.to_csv(keep_tsv, sep="\t", index=False)

    pd.DataFrame([{
        "sample": os.path.basename(args.reditools_bed),
        "min_cov_q30": args.min_cov,
        "min_edited_sites_per_interval": args.min_edited_sites_per_interval,
        "n_intervals_total_with_edits": int(outdf.shape[0]),
        "n_intervals_used": int(keep.shape[0]),
        "pearson_r": pear,
        "spearman_r": spear,
        "n_points_corr": n,
        "y_used": "AG_edit_fraction"
    }]).to_csv(stats_tsv, sep="\t", index=False)

    # -----------------------------
    # Plot 1: binned trend curve
    # -----------------------------
    binned_tsv = args.out_prefix + ".binned_curve.tsv"
    curve_png = args.out_prefix + ".binned_curve.png"

    if keep.shape[0] > 0:
        nbins = args.bins if args.bins and args.bins > 0 else (30 if keep.shape[0] >= 300 else 15)
        keep2 = keep.copy()
        keep2["force_bin"] = safe_qcut(keep2["max_mean_dsRNA_force"], q=nbins)

        b = (keep2
             .groupby("force_bin", observed=False)
             .agg(
                 force_mid=("max_mean_dsRNA_force", "median"),
                 ag_frac_mean=("AG_edit_fraction", "mean"),
                 ag_frac_median=("AG_edit_fraction", "median"),
                 n=("AG_edit_fraction", "size"),
             )
             .reset_index()
             .sort_values("force_mid"))

        b.to_csv(binned_tsv, sep="\t", index=False)

        plt.figure()
        plt.plot(b["force_mid"], b["ag_frac_mean"], marker="o")
        plt.xlabel("max_mean_dsRNA_force (repeat interval)")
        plt.ylabel("Fraction of edited sites that are AG")
        plt.title(
            f"{os.path.basename(args.reditools_bed)}\n"
            f"Binned curve (nbins={nbins}) | Pearson r={pear:.3g} Spearman ρ={spear:.3g} n={n}"
        )
        plt.tight_layout()
        plt.savefig(curve_png, dpi=200)
        plt.close()
    else:
        pd.DataFrame(columns=["force_mid", "ag_frac_mean", "ag_frac_median", "n"]).to_csv(binned_tsv, sep="\t", index=False)

    # -----------------------------
    # Plot 2: density curves by force tertile
    # -----------------------------
    dens_png = args.out_prefix + ".density_by_force_group.png"
    if keep.shape[0] > 0:
        keep3 = keep.copy()
        keep3["force_group"] = safe_qcut(keep3["max_mean_dsRNA_force"], q=3, labels=["low", "mid", "high"])

        plt.figure()
        for grp in ["low", "mid", "high"]:
            vals = keep3.loc[keep3["force_group"] == grp, "AG_edit_fraction"].dropna().values
            if len(vals) < 5:
                continue
            hist, edges = np.histogram(vals, bins=50, density=True)
            mids = (edges[:-1] + edges[1:]) / 2
            plt.plot(mids, hist, label=f"{grp} force (n={len(vals)})")

        plt.xlabel("AG_edit_fraction (per interval)")
        plt.ylabel("Density")
        plt.title(f"{os.path.basename(args.reditools_bed)}\nAG site-fraction distribution by force group")
        plt.legend()
        plt.tight_layout()
        plt.savefig(dens_png, dpi=200)
        plt.close()

    print(
        "[OK] Wrote:\n"
        f"  {out_tsv}\n"
        f"  {keep_tsv}\n"
        f"  {stats_tsv}\n"
        f"  {binned_tsv}\n"
        f"  {curve_png}\n"
        f"  {dens_png}\n"
    )


if __name__ == "__main__":
    main()

