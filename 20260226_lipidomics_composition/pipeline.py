"""
Lipidomics Composition Analysis Pipeline
Usage:
    python pipeline.py SM
    python pipeline.py PC
    python pipeline.py PE
    python pipeline.py PI
    python pipeline.py PS
    python pipeline.py Cer
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# ---------- 設定 ----------
TREATMENTS   = ["Ctr", "C18:0"]
GROUPS       = ["Resistant", "Sensitive"]
DATA_DIR     = "data"
RESULTS_DIR  = "results"


def load_data(lipid_class):
    path = os.path.join(DATA_DIR, f"{lipid_class}_pct_individual.csv")
    if not os.path.exists(path):
        print(f"[ERROR] データファイルが見つかりません: {path}")
        sys.exit(1)
    df = pd.read_csv(path)
    print(f"[{lipid_class}] データ読み込み完了: {df.shape[0]} 行, species {df['species'].nunique()} 種")
    return df


def out_dir(lipid_class):
    d = os.path.join(RESULTS_DIR, lipid_class)
    os.makedirs(d, exist_ok=True)
    return d


# ---- Figure 1: Stacked bar (Resistant vs Sensitive) ----
def plot_stacked_bar(df, lipid_class, treatment, outdir):
    sub = df[df["treatment"] == treatment]
    summary = sub.groupby(["group", "species"])["pct"].mean().reset_index()
    pivot = summary.pivot(index="species", columns="group", values="pct")
    pivot = pivot.sort_values("Resistant", ascending=False)

    n = len(pivot)
    colors = cm.tab20(np.linspace(0, 1, n))

    fig, ax = plt.subplots(figsize=(6, 6))
    bottoms = np.zeros(2)
    for i, sp in enumerate(pivot.index):
        vals = [pivot.loc[sp, g] for g in GROUPS]
        ax.bar(GROUPS, vals, bottom=bottoms, color=colors[i], label=sp, width=0.5)
        bottoms += np.array(vals)

    ax.set_ylabel("mol% of total " + lipid_class, fontsize=12)
    ax.set_title(f"{lipid_class} Species Composition\n(Resistant vs Sensitive, {treatment})", fontsize=12)
    ax.set_ylim(0, 105)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8,
              frameon=False, title="Species", title_fontsize=9)
    plt.tight_layout()

    tag = treatment.replace(":", "").replace("+", "_")
    stem = os.path.join(outdir, f"{lipid_class}_composition_{tag}")
    plt.savefig(stem + ".pdf", bbox_inches="tight")
    plt.savefig(stem + ".png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {stem}.png")


# ---- Figure 2: Volcano - Resistant vs Sensitive ----
def plot_volcano_res_vs_sen(df, lipid_class, outdir):
    species_list = sorted(df["species"].unique())
    records = []
    for trt in TREATMENTS:
        sub = df[df["treatment"] == trt]
        for sp in species_list:
            res = sub[(sub["group"] == "Resistant") & (sub["species"] == sp)]["pct"].values
            sen = sub[(sub["group"] == "Sensitive") & (sub["species"] == sp)]["pct"].values
            if len(res) < 3 or len(sen) < 3:
                continue
            _, p = mannwhitneyu(res, sen, alternative="two-sided")
            fc = np.log2((res.mean() + 1e-6) / (sen.mean() + 1e-6))
            records.append({"treatment": trt, "species": sp, "p": p, "log2FC": fc})

    result = pd.DataFrame(records)
    for trt in TREATMENTS:
        idx = result["treatment"] == trt
        if idx.sum() > 1:
            _, q, _, _ = multipletests(result.loc[idx, "p"], method="fdr_bh")
            result.loc[idx, "q"] = q
    result["-log10p"] = -np.log10(result["p"] + 1e-10)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, trt in zip(axes, TREATMENTS):
        sub = result[result["treatment"] == trt].copy()

        colors = sub.apply(lambda r: "#e05a5a" if r["p"] < 0.05 and r["log2FC"] > 0
                           else ("#4a90d9" if r["p"] < 0.05 and r["log2FC"] < 0
                                 else "#cccccc"), axis=1)
        ax.scatter(sub["log2FC"], sub["-log10p"], c=colors, s=55, alpha=0.85, edgecolors="none")
        ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8)

        for _, row in sub[sub["p"] < 0.05].iterrows():
            ax.text(row["log2FC"], row["-log10p"] + 0.04, row["species"],
                    fontsize=6.5, ha="center", va="bottom")

        sig = sub[sub.get("q", pd.Series([1]*len(sub))) < 0.2]["species"].tolist()
        note = "FDR q<0.2: " + ", ".join(sig) if sig else "no significance after FDR"
        ax.set_xlabel("log2 FC (Resistant / Sensitive)", fontsize=11)
        ax.set_ylabel("-log10(p-value)", fontsize=11)
        ax.set_title(f"{lipid_class}: Resistant vs Sensitive — {trt}\n({note})", fontsize=10)
        patches = [mpatches.Patch(color="#e05a5a", label="Higher in Resistant (p<0.05)"),
                   mpatches.Patch(color="#4a90d9", label="Higher in Sensitive (p<0.05)"),
                   mpatches.Patch(color="#cccccc", label="n.s.")]
        ax.legend(handles=patches, fontsize=8, frameon=False)

    plt.tight_layout()
    stem = os.path.join(outdir, f"{lipid_class}_volcano_ResVsSen")
    plt.savefig(stem + ".pdf", bbox_inches="tight")
    plt.savefig(stem + ".png", dpi=150, bbox_inches="tight")
    plt.close()

    result.to_csv(os.path.join(outdir, f"{lipid_class}_stats_ResVsSen.csv"), index=False)
    print(f"  Saved: {stem}.png")


# ---- Figure 3: Volcano - Ctr vs C18:0 (treatment effect) ----
def plot_volcano_treatment(df, lipid_class, outdir):
    species_list = sorted(df["species"].unique())
    records = []
    for grp in GROUPS:
        sub = df[df["group"] == grp]
        for sp in species_list:
            ctr = sub[(sub["treatment"] == "Ctr")   & (sub["species"] == sp)]["pct"].values
            c18 = sub[(sub["treatment"] == "C18:0") & (sub["species"] == sp)]["pct"].values
            if len(ctr) < 3 or len(c18) < 3:
                continue
            _, p = mannwhitneyu(ctr, c18, alternative="two-sided")
            fc = np.log2((c18.mean() + 1e-6) / (ctr.mean() + 1e-6))
            records.append({"group": grp, "species": sp, "p": p, "log2FC": fc})

    result = pd.DataFrame(records)
    for grp in GROUPS:
        idx = result["group"] == grp
        if idx.sum() > 1:
            _, q, _, _ = multipletests(result.loc[idx, "p"], method="fdr_bh")
            result.loc[idx, "q"] = q
    result["-log10p"] = -np.log10(result["p"] + 1e-10)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, grp in zip(axes, GROUPS):
        sub = result[result["group"] == grp].copy()

        colors = sub.apply(lambda r: "#e07030" if r["p"] < 0.05 and r["log2FC"] > 0
                           else ("#5090d0" if r["p"] < 0.05 and r["log2FC"] < 0
                                 else "#cccccc"), axis=1)
        ax.scatter(sub["log2FC"], sub["-log10p"], c=colors, s=55, alpha=0.85, edgecolors="none")
        ax.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8)

        for _, row in sub[sub["p"] < 0.05].iterrows():
            ax.text(row["log2FC"], row["-log10p"] + 0.04, row["species"],
                    fontsize=6.5, ha="center", va="bottom")

        sig = sub[sub.get("q", pd.Series([1]*len(sub))) < 0.2]["species"].tolist()
        note = "FDR q<0.2: " + ", ".join(sig) if sig else "no significance after FDR"
        ax.set_xlabel("log2 FC (C18:0 / Ctr)", fontsize=11)
        ax.set_ylabel("-log10(p-value)", fontsize=11)
        ax.set_title(f"{lipid_class}: {grp} — Ctr vs C18:0\n({note})", fontsize=10)
        patches = [mpatches.Patch(color="#e07030", label="Higher in C18:0 (p<0.05)"),
                   mpatches.Patch(color="#5090d0", label="Lower in C18:0 (p<0.05)"),
                   mpatches.Patch(color="#cccccc", label="n.s.")]
        ax.legend(handles=patches, fontsize=8, frameon=False)

    plt.tight_layout()
    stem = os.path.join(outdir, f"{lipid_class}_volcano_treatment")
    plt.savefig(stem + ".pdf", bbox_inches="tight")
    plt.savefig(stem + ".png", dpi=150, bbox_inches="tight")
    plt.close()

    result.to_csv(os.path.join(outdir, f"{lipid_class}_stats_treatment.csv"), index=False)
    print(f"  Saved: {stem}.png")


# ---------- メイン ----------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python pipeline.py <LipidClass>")
        print("Example: python pipeline.py SM")
        sys.exit(1)

    lipid_class = sys.argv[1]
    print(f"\n===== {lipid_class} 解析開始 =====")

    df = load_data(lipid_class)
    outdir = out_dir(lipid_class)

    print("[1/3] Stacked bar charts ...")
    for trt in TREATMENTS:
        plot_stacked_bar(df, lipid_class, trt, outdir)

    print("[2/3] Volcano: Resistant vs Sensitive ...")
    plot_volcano_res_vs_sen(df, lipid_class, outdir)

    print("[3/3] Volcano: treatment effect (Ctr vs C18:0) ...")
    plot_volcano_treatment(df, lipid_class, outdir)

    print(f"\n===== 完了: 出力先 {outdir}/ =====\n")
