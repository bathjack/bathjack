"""
Lipidomics Composition Analysis Pipeline
Usage:
    python pipeline.py SM
    python pipeline.py SM --experiment isotope
    python pipeline.py PC --experiment unlabeled
"""

import sys
import os
import argparse
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
TREATMENTS   = ["Ctr", "C18:0", "C18:0+C18:1"]
GROUPS       = ["Resistant", "Sensitive"]
BASE_DATA    = "data"
BASE_RESULTS = "results"


def load_data(lipid_class, experiment):
    path = os.path.join(BASE_DATA, experiment, f"{lipid_class}_pct_individual.csv")
    if not os.path.exists(path):
        print(f"[ERROR] データファイルが見つかりません: {path}")
        sys.exit(1)
    df = pd.read_csv(path)
    print(f"[{lipid_class}] データ読み込み完了: {df.shape[0]} 行, species {df['species'].nunique()} 種")
    return df


def out_dir(lipid_class, experiment):
    d = os.path.join(BASE_RESULTS, experiment, lipid_class)
    os.makedirs(d, exist_ok=True)
    return d


# ---- Figure 1: Stacked bar (Resistant vs Sensitive) ----
OTHERS_THRESHOLD = 1.0   # 全グループの平均 mol% がこの値未満の種を "Others" にまとめる

def plot_stacked_bar(df, lipid_class, treatment, outdir):
    sub = df[df["treatment"] == treatment]
    summary = sub.groupby(["group", "species"])["pct"].mean().reset_index()
    pivot = summary.pivot(index="species", columns="group", values="pct").fillna(0)
    pivot = pivot.sort_values("Resistant", ascending=False)

    # 閾値未満の種を Others にまとめる
    max_pct = pivot.max(axis=1)
    major = pivot[max_pct >= OTHERS_THRESHOLD]
    minor = pivot[max_pct <  OTHERS_THRESHOLD]
    if len(minor) > 0:
        others_row = minor.sum(axis=0).rename("Others")
        pivot_plot = pd.concat([major, others_row.to_frame().T])
    else:
        pivot_plot = major

    n = len(pivot_plot)
    # Others はグレー、それ以外は tab20
    colors = list(cm.tab20(np.linspace(0, 1, n - (1 if len(minor) > 0 else 0))))
    if len(minor) > 0:
        colors.append((0.75, 0.75, 0.75, 1.0))   # gray for Others

    fig, ax = plt.subplots(figsize=(6, 6))
    bottoms = np.zeros(2)
    for i, sp in enumerate(pivot_plot.index):
        vals = [pivot_plot.loc[sp, g] for g in GROUPS]
        ax.bar(GROUPS, vals, bottom=bottoms, color=colors[i], label=sp, width=0.5)
        bottoms += np.array(vals)

    ax.set_ylabel("mol% of total " + lipid_class, fontsize=12)
    ax.set_title(f"{lipid_class} Species Composition\n(Resistant vs Sensitive, {treatment})", fontsize=12)
    ax.set_ylim(0, 105)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8,
              frameon=False, title=f"Species (≥{OTHERS_THRESHOLD}%)", title_fontsize=9)
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

        label_df = sub[sub["p"] < 0.05].nsmallest(10, "p")
        for _, row in label_df.iterrows():
            ax.text(row["log2FC"], row["-log10p"] + 0.04, row["species"],
                    fontsize=6.5, ha="center", va="bottom")

        sig = sub[sub.get("q", pd.Series([1]*len(sub))) < 0.2]["species"].tolist()
        note = f"FDR q<0.2: {len(sig)} species" if sig else "no significance after FDR"
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

        label_df = sub[sub["p"] < 0.05].nsmallest(10, "p")
        for _, row in label_df.iterrows():
            ax.text(row["log2FC"], row["-log10p"] + 0.04, row["species"],
                    fontsize=6.5, ha="center", va="bottom")

        sig = sub[sub.get("q", pd.Series([1]*len(sub))) < 0.2]["species"].tolist()
        note = f"FDR q<0.2: {len(sig)} species" if sig else "no significance after FDR"
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


# ---- Figure 4: Rescue scatter (C18:0 effect vs C18:0+C18:1 effect) ----
def plot_rescue(df, lipid_class, outdir):
    """
    C18:1 の rescue 効果を可視化するスキャッタープロット。
      x軸: log2FC (C18:0 / Ctr)          ← C18:0 の変化
      y軸: log2FC (C18:0+C18:1 / Ctr)    ← C18:0+C18:1 の変化
    対角線 (y=x) より y が小さい → C18:1 が C18:0 の効果を打ち消している
    """
    RESCUE_TRT = "C18:0+C18:1"
    if RESCUE_TRT not in df["treatment"].unique():
        print(f"  [{lipid_class}] '{RESCUE_TRT}' データなし — スキップ")
        return

    species_list = sorted(df["species"].unique())
    records = []
    for grp in GROUPS:
        sub = df[df["group"] == grp]
        for sp in species_list:
            ctr = sub[(sub["treatment"] == "Ctr")      & (sub["species"] == sp)]["pct"].values
            c18 = sub[(sub["treatment"] == "C18:0")    & (sub["species"] == sp)]["pct"].values
            rsc = sub[(sub["treatment"] == RESCUE_TRT) & (sub["species"] == sp)]["pct"].values
            if len(ctr) < 3 or len(c18) < 3 or len(rsc) < 3:
                continue
            _, p_c18 = mannwhitneyu(ctr, c18, alternative="two-sided")
            _, p_rsc = mannwhitneyu(c18, rsc, alternative="two-sided")
            fc_c18 = np.log2((c18.mean() + 1e-6) / (ctr.mean() + 1e-6))
            fc_rsc = np.log2((rsc.mean() + 1e-6) / (ctr.mean() + 1e-6))
            records.append({
                "group": grp, "species": sp,
                "log2FC_C18": fc_c18, "log2FC_rescue": fc_rsc,
                "p_C18": p_c18, "p_C18vsRescue": p_rsc,
                "mean_Ctr": ctr.mean(), "mean_C18": c18.mean(), "mean_rescue": rsc.mean(),
            })

    if not records:
        print(f"  [{lipid_class}] rescue データ不足 — スキップ")
        return

    result = pd.DataFrame(records)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, grp in zip(axes, GROUPS):
        sub = result[result["group"] == grp].copy()
        if len(sub) == 0:
            continue

        # C18:0 で有意に変化したものを色付け
        colors = sub.apply(
            lambda r: "#e07030" if r["p_C18"] < 0.05 and r["log2FC_C18"] > 0
                      else ("#5090d0" if r["p_C18"] < 0.05 and r["log2FC_C18"] < 0
                            else "#cccccc"), axis=1
        )
        ax.scatter(sub["log2FC_C18"], sub["log2FC_rescue"],
                   c=colors, s=55, alpha=0.85, edgecolors="none")

        # 対角線 y=x (no rescue)
        lim = max(sub["log2FC_C18"].abs().max(), sub["log2FC_rescue"].abs().max()) * 1.15
        ax.plot([-lim, lim], [-lim, lim], "k--", linewidth=0.9, alpha=0.5, label="y = x (no rescue)")
        ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.8)

        # C18:0 で有意変化した上位 10 種にラベル
        label_df = sub[sub["p_C18"] < 0.05].nsmallest(10, "p_C18")
        for _, row in label_df.iterrows():
            ax.text(row["log2FC_C18"], row["log2FC_rescue"] + 0.04, row["species"],
                    fontsize=6.5, ha="center", va="bottom")

        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("log2 FC (C18:0 / Ctr)", fontsize=11)
        ax.set_ylabel("log2 FC (C18:0+C18:1 / Ctr)", fontsize=11)
        ax.set_title(f"{lipid_class}: {grp}\nRescue by C18:1", fontsize=11)

        patches = [mpatches.Patch(color="#e07030", label="Increased by C18:0 (p<0.05)"),
                   mpatches.Patch(color="#5090d0", label="Decreased by C18:0 (p<0.05)"),
                   mpatches.Patch(color="#cccccc", label="n.s.")]
        ax.legend(handles=patches, fontsize=8, frameon=False)
        ax.text(0.02, 0.98, "← rescued", transform=ax.transAxes,
                fontsize=8, color="#5090d0", va="top", style="italic")
        ax.text(0.98, 0.02, "exacerbated →", transform=ax.transAxes,
                fontsize=8, color="#e07030", va="bottom", ha="right", style="italic")

    plt.tight_layout()
    stem = os.path.join(outdir, f"{lipid_class}_rescue_C18")
    plt.savefig(stem + ".pdf", bbox_inches="tight")
    plt.savefig(stem + ".png", dpi=150, bbox_inches="tight")
    plt.close()

    result.to_csv(os.path.join(outdir, f"{lipid_class}_stats_rescue.csv"), index=False)
    print(f"  Saved: {stem}.png")


# ---------- メイン ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lipidomics composition pipeline")
    parser.add_argument("lipid_class",
                        help="脂質クラス名（SM, Cer, PC, PE, PI, PS）")
    parser.add_argument("--experiment", default="unlabeled",
                        help="実験名サブディレクトリ（デフォルト: unlabeled）")
    args = parser.parse_args()

    lipid_class = args.lipid_class
    experiment  = args.experiment

    print(f"\n===== {lipid_class} 解析開始 [{experiment}] =====")

    df = load_data(lipid_class, experiment)
    outdir = out_dir(lipid_class, experiment)

    avail_trts = df["treatment"].unique()

    print("[1/3] Stacked bar charts ...")
    for trt in TREATMENTS:
        if trt in avail_trts:
            plot_stacked_bar(df, lipid_class, trt, outdir)

    print("[2/3] Volcano: Resistant vs Sensitive ...")
    plot_volcano_res_vs_sen(df, lipid_class, outdir)

    print("[3/3] Volcano: treatment effect (Ctr vs C18:0) ...")
    plot_volcano_treatment(df, lipid_class, outdir)

    print("[4/4] Rescue plot: C18:1 protection ...")
    plot_rescue(df, lipid_class, outdir)

    print(f"\n===== 完了: 出力先 {outdir}/ =====\n")
