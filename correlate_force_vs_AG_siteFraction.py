#!/usr/bin/env python3
import os
import re
import sys
import argparse
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from intervaltree import IntervalTree


# -----------------------------
# Helpers
# -----------------------------
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def read_tsv_auto(path: str) -> pd.DataFrame:
    # Robust TSV/CSV reader
    return pd.read_csv(path, sep=None, engine="python", low_memory=True)


def build_interval_trees_from_bed(bed_path: str):
    """
    BED: chr start end [name ...]
    """
    trees = {}
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            if chrom not in trees:
                trees[chrom] = IntervalTree()
            trees[chrom].addi(start, end, True)
    return trees


def overlaps_any(trees, chrom, pos0):
    """
    pos0: 0-based coordinate. checks if pos0 in any interval [start,end)
    """
    if chrom not in trees:
        return False
    return bool(trees[chrom].overlap(pos0, pos0 + 1))


def build_force_intervaltrees(force_tsv: str):
    """
    Expected columns include:
      chr, start, end, max_mean_dsRNA_force, class, family
    Uses *lowercase* class/family.
    """
    df = pd.read_csv(force_tsv, sep="\t")

    required = {"chr", "start", "end", "max_mean_dsRNA_force", "class", "family"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"[FATAL] Missing required columns in force TSV: {sorted(list(missing))}\n"
            f"Columns found: {list(df.columns)}"
        )

    df["start"] = pd.to_numeric(df["start"])
    df["end"] = pd.to_numeric(df["end"])
    df["max_mean_dsRNA_force"] = pd.to_numeric(df["max_mean_dsRNA_force"], errors="coerce")

    trees = {}
    meta = {}

    for _, r in df.iterrows():
        chrom = str(r["chr"])
        s = int(r["start"])
        e = int(r["end"])

        rep_id = str(r.get("rep.id", "NA"))
        strand = str(r.get("strand", "."))
        key = f"{chrom}:{s}-{e}|{rep_id}|{strand}"

        if chrom not in trees:
            trees[chrom] = IntervalTree()
        trees[chrom].addi(s, e, key)

        meta[key] = {
            "chr": chrom,
            "start": s,
            "end": e,
            "rep.id": rep_id,
            "strand": strand,
            "class": str(r["class"]),      # LOWERCASE
            "family": str(r["family"]),    # LOWERCASE
            "max_mean_dsRNA_force": float(r["max_mean_dsRNA_force"]) if pd.notnull(r["max_mean_dsRNA_force"]) else np.nan,
        }

    return trees, meta


def infer_reditools_columns(df: pd.DataFrame):
    """
    Tries to find:
      chrom column, pos column, substitution column, editing-index column
    Your example looks like it already has:
      Region/chr, Position, AllSubs (AG), and a frequency like 0.22
    """

    colmap = {c.lower(): c for c in df.columns}

    # Chrom
    chrom_candidates = ["region", "chr", "chrom", "chromosome"]
    chrom_col = next((colmap[c] for c in chrom_candidates if c in colmap), None)

    # Position (1-based)
    pos_candidates = ["position", "pos", "start"]
    pos_col = next((colmap[c] for c in pos_candidates if c in colmap), None)

    # Substitution type
    allsubs_candidates = ["allsubs", "substitution", "subs"]
    allsubs_col = next((colmap[c] for c in allsubs_candidates if c in colmap), None)

    # Editing index / frequency
    freq_candidates = ["frequency", "freq", "editing_index", "editindex"]
    freq_col = next((colmap[c] for c in freq_candidates if c in colmap), None)

    if chrom_col is None or pos_col is None:
        raise ValueError(
            f"[FATAL] Could not infer chrom/pos columns from: {list(df.columns)}.\n"
            f"Need something like 'Region/chr' and 'Position'."
        )

    if allsubs_col is None:
        # Sometimes it's literally 'AllSubs' but not found; user file might have 'AllSubs' with capitalization
        # If missing entirely, we can still proceed but only if there's an AG-specific column.
        eprint("[WARN] Could not find AllSubs/substitution column; will try to proceed without it.")
    if freq_col is None:
        eprint("[WARN] Could not find Frequency/Freq column; will try to compute editing index from AG counts if available.")

    return chrom_col, pos_col, allsubs_col, freq_col


def compute_editing_index(df: pd.DataFrame, allsubs_col: str | None, freq_col: str | None):
    """
    Returns:
      chrom (str), pos1 (int), is_AG (bool), editing_index (float)
    If freq_col exists, uses that as editing_index.
    Otherwise tries to compute from counts if possible (very file-dependent).
    """
    out = pd.DataFrame()

    chrom_col, pos_col, allsubs_col2, freq_col2 = infer_reditools_columns(df)
    allsubs_col = allsubs_col or allsubs_col2
    freq_col = freq_col or freq_col2

    out["chr"] = df[chrom_col].astype(str)
    out["pos1"] = pd.to_numeric(df[pos_col], errors="coerce").astype("Int64")

    if allsubs_col and allsubs_col in df.columns:
        out["allsubs"] = df[allsubs_col].astype(str)
        out["is_AG"] = out["allsubs"].str.upper().str.contains(r"\bAG\b") | (out["allsubs"].str.upper() == "AG")
    else:
        out["allsubs"] = ""
        out["is_AG"] = False

    if freq_col and freq_col in df.columns:
        out["editing_index"] = pd.to_numeric(df[freq_col], errors="coerce")
    else:
        # Fallback heuristic: look for a column that seems like AG frequency
        guess = None
        for c in df.columns:
            cl = c.lower()
            if "ag" in cl and ("freq" in cl or "frequency" in cl or "rate" in cl):
                guess = c
                break
        if guess:
            out["editing_index"] = pd.to_numeric(df[guess], errors="coerce")
        else:
            # last-resort: NaN
            out["editing_index"] = np.nan

    out = out.dropna(subset=["pos1"])
    out["pos0"] = out["pos1"].astype(int) - 1  # convert to 0-based
    return out


# -----------------------------
# Plotting
# -----------------------------
def violin_plot(groups, out_png, title):
    """
    groups: dict name -> array-like values
    """
    names = [k for k, v in groups.items() if len(v) > 0]
    data = [np.asarray(groups[k], dtype=float) for k in names]

    plt.figure(figsize=(8, 5))
    plt.violinplot(data, showmeans=True, showextrema=True, showmedians=True)
    plt.xticks(range(1, len(names) + 1), names, rotation=20, ha="right")
    plt.ylabel("A→G editing index")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def density_plot(force_all, force_high, out_png, title, high_label="high-editing"):
    plt.figure(figsize=(8, 5))
    # Use histogram density (robust on clusters)
    bins = 60
    force_all = np.asarray(force_all, dtype=float)
    force_high = np.asarray(force_high, dtype=float)

    force_all = force_all[np.isfinite(force_all)]
    force_high = force_high[np.isfinite(force_high)]

    if len(force_all) == 0:
        eprint("[WARN] No force values to plot density.")
        return

    plt.hist(force_all, bins=bins, density=True, alpha=0.5, label="all mapped sites")
    if len(force_high) > 0:
        plt.hist(force_high, bins=bins, density=True, alpha=0.5, label=high_label)

    plt.xlabel("max_mean_dsRNA_force")
    plt.ylabel("Density")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def scatter_plot(df, out_png, title):
    """
    df must have: max_mean_dsRNA_force, editing_index, group
    """
    d = df.copy()
    d = d[np.isfinite(d["max_mean_dsRNA_force"]) & np.isfinite(d["editing_index"])]

    if len(d) == 0:
        eprint("[WARN] No points to plot scatter.")
        return

    # Downsample if huge
    if len(d) > 200000:
        d = d.sample(n=200000, random_state=1)

    plt.figure(figsize=(7, 6))
    # plot each group separately (no custom colors requested; matplotlib will assign)
    for gname, sub in d.groupby("group", dropna=False):
        plt.scatter(sub["max_mean_dsRNA_force"], sub["editing_index"], s=6, alpha=0.3, label=str(gname))

    plt.xlabel("max_mean_dsRNA_force (dsRNA force)")
    plt.ylabel("A→G editing index")
    plt.title(title)
    plt.legend(markerscale=2, fontsize=9)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Intersect REDItools A→G sites with dsRNA force intervals; plot MOR vs IR/nonIR Alu violins and force vs editing."
    )
    ap.add_argument("--reditools", required=True, help="REDItools TSV for ONE sample")
    ap.add_argument("--force_tsv", required=True, help="squire_seqA_intersect_100%_annot.tsv")
    ap.add_argument("--bed_mor", required=True, help="MOR BED (chr start end)")
    ap.add_argument("--bed_ir_alu", required=True, help="IR Alu BED (chr start end)")
    ap.add_argument("--bed_nonir_alu", required=True, help="non-IR Alu BED (chr start end)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--sample", default=None, help="Sample name override (default = basename of reditools file)")
    ap.add_argument("--high_edit_thresh", type=float, default=0.10, help="Editing index threshold for 'high editing' overlay (default 0.10)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    sample = args.sample
    if not sample:
        sample = os.path.basename(args.reditools)
        sample = re.sub(r"\.(tsv|txt|csv)(\.gz)?$", "", sample)

    eprint(f"[INFO] Sample: {sample}")
    eprint(f"[INFO] Reading REDItools: {args.reditools}")
    redi_raw = read_tsv_auto(args.reditools)

    eprint(f"[INFO] Reading force TSV: {args.force_tsv}")
    force_trees, force_meta = build_force_intervaltrees(args.force_tsv)

    eprint("[INFO] Reading BEDs (MOR, IR_Alu, nonIR_Alu)")
    mor_trees = build_interval_trees_from_bed(args.bed_mor)
    ir_trees = build_interval_trees_from_bed(args.bed_ir_alu)
    nonir_trees = build_interval_trees_from_bed(args.bed_nonir_alu)

    # Compute editing index / A→G flags
    redi = compute_editing_index(redi_raw, allsubs_col=None, freq_col=None)
    redi = redi.dropna(subset=["editing_index"])
    eprint(f"[INFO] Sites with editing_index: {len(redi):,}")

    # Map each site to force interval (first overlap; if multiple, choose max force)
    mapped_rows = []
    n_mapped = 0

    for r in redi.itertuples(index=False):
        chrom = r.chr
        pos0 = int(r.pos0)

        if chrom not in force_trees:
            continue

        hits = force_trees[chrom].overlap(pos0, pos0 + 1)
        if not hits:
            continue

        # choose interval with max dsRNA force (more robust than arbitrary first)
        best_key = None
        best_force = -np.inf
        for h in hits:
            key = h.data
            f = force_meta[key]["max_mean_dsRNA_force"]
            if np.isfinite(f) and f > best_force:
                best_force = f
                best_key = key

        if best_key is None:
            # all NaN forces; still keep first hit
            best_key = next(iter(hits)).data
            best_force = force_meta[best_key]["max_mean_dsRNA_force"]

        m = force_meta[best_key]

        # Group assignment by BED overlap (priority MOR > IR > nonIR)
        group = None
        if overlaps_any(mor_trees, chrom, pos0):
            group = "MOR"
        elif overlaps_any(ir_trees, chrom, pos0):
            group = "IR_Alu"
        elif overlaps_any(nonir_trees, chrom, pos0):
            group = "nonIR_Alu"
        else:
            continue

        n_mapped += 1

        mapped_rows.append({
            "sample": sample,
            "chr": chrom,
            "pos1": int(r.pos1),
            "pos0": pos0,
            "is_AG": bool(r.is_AG),
            "editing_index": float(r.editing_index),
            "group": group,
            # force interval info
            "force_chr": m["chr"],
            "force_start": m["start"],
            "force_end": m["end"],
            "rep.id": m["rep.id"],
            "strand": m["strand"],
            "class": m["class"],     # lowercase
            "family": m["family"],   # lowercase
            "max_mean_dsRNA_force": m["max_mean_dsRNA_force"],
        })

    eprint(f"[INFO] Mapped+grouped sites: {n_mapped:,}")

    if n_mapped == 0:
        raise SystemExit("[FATAL] No sites mapped to force intervals AND assigned to MOR/IR/nonIR. Check your BEDs + chromosome naming (chr1 vs 1).")

    merged = pd.DataFrame(mapped_rows)

    # Save merged table
    out_tsv = os.path.join(args.outdir, f"{sample}.merged_force_editing.tsv")
    merged.to_csv(out_tsv, sep="\t", index=False)
    eprint(f"[INFO] Wrote merged TSV: {out_tsv}")

    # -----------------
    # Plot 1: Violin (editing index by group)
    # -----------------
    groups = {
        "MOR": merged.loc[merged["group"] == "MOR", "editing_index"].dropna().values,
        "IR_Alu": merged.loc[merged["group"] == "IR_Alu", "editing_index"].dropna().values,
        "nonIR_Alu": merged.loc[merged["group"] == "nonIR_Alu", "editing_index"].dropna().values,
    }
    out_violin = os.path.join(args.outdir, f"{sample}.violin_editing_index.png")
    violin_plot(groups, out_violin, f"{sample}: A→G editing index (MOR vs IR/nonIR Alu)")
    eprint(f"[INFO] Wrote: {out_violin}")

    # -----------------
    # Plot 2: Force density (all mapped vs high-editing)
    # -----------------
    force_all = merged["max_mean_dsRNA_force"].values
    force_high = merged.loc[merged["editing_index"] >= args.high_edit_thresh, "max_mean_dsRNA_force"].values
    out_dens = os.path.join(args.outdir, f"{sample}.force_density.png")
    density_plot(
        force_all, force_high, out_dens,
        f"{sample}: dsRNA force density (all vs editing≥{args.high_edit_thresh})",
        high_label=f"editing≥{args.high_edit_thresh}"
    )
    eprint(f"[INFO] Wrote: {out_dens}")

    # -----------------
    # Plot 3: Scatter force vs editing (colored by group)
    # -----------------
    out_scatter = os.path.join(args.outdir, f"{sample}.force_vs_editing_scatter.png")
    scatter_plot(merged, out_scatter, f"{sample}: dsRNA force vs A→G editing index")
    eprint(f"[INFO] Wrote: {out_scatter}")

    # quick stats
    stats_out = os.path.join(args.outdir, f"{sample}.summary.txt")
    with open(stats_out, "w") as w:
        w.write(f"sample\t{sample}\n")
        w.write(f"n_mapped_grouped\t{len(merged)}\n")
        for g in ["MOR", "IR_Alu", "nonIR_Alu"]:
            sub = merged[merged["group"] == g]
            w.write(f"{g}_n\t{len(sub)}\n")
            w.write(f"{g}_median_editing\t{sub['editing_index'].median() if len(sub) else np.nan}\n")
            w.write(f"{g}_median_force\t{sub['max_mean_dsRNA_force'].median() if len(sub) else np.nan}\n")
    eprint(f"[INFO] Wrote: {stats_out}")


if __name__ == "__main__":
    main()
