"""
Lipidomics Isotope Experiment Pipeline
Data: isotopeFFA_vitro_summary.xlsx (tidy_percent sheet)

Cell line grouping:
  Sensitive (lipotoxicity+): OVCAR5, OVCAR8, SKOV3, ES2, OVCAR3 (ovarian cancer), HOSE (normal ovarian)
  Resistant (lipotoxicity-): H1299 (lung cancer), MCF10A (normal mammary epithelial)

Usage:
    python pipeline_isotope.py PI
    python pipeline_isotope.py PC
    python pipeline_isotope.py PI --experiment isotope
    python pipeline_isotope.py PI --focus "PI(18:0_18:0)"
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

# ---------- 設定 ----------
SENSITIVE = ["OVCAR5", "OVCAR8", "SKOV3", "ES2", "OVCAR3", "HOSE"]
RESISTANT = ["H1299", "MCF10A"]
GROUPS    = {"Sensitive": SENSITIVE, "Resistant": RESISTANT}
GROUP_COLORS = {"Sensitive": "#e05a5a", "Resistant": "#4a90d9"}

TREATMENTS = ["BSA", "18:0-d35", "18:1-13C5", "18:0-d35+18:1-13C5"]
TREATMENT_TAGS = {
    "BSA":                  "BSA",
    "18:0-d35":             "18_0-d35",
    "18:1-13C5":            "18_1-13C5",
    "18:0-d35+18:1-13C5":  "combo",
}

BASE_DATA     = "data"
BASE_RESULTS  = "results"
EXCEL_FILE    = "isotopeFFA_vitro_summary.xlsx"
DEFAULT_SHEET = "tidy_percent"
OTHERS_THRESHOLD = 1.0


def load_data(lipid_class, experiment="isotope", sheet=DEFAULT_SHEET):
    path = os.path.join(BASE_DATA, experiment, EXCEL_FILE)
    if not os.path.exists(path):
        print(f"[ERROR] ファイルが見つかりません: {path}")
        sys.exit(1)

    df = pd.read_excel(path, sheet_name=sheet)
    df = df[df["class"] == lipid_class].copy()

    if len(df) == 0:
        print(f"[ERROR] クラス '{lipid_class}' のデータが見つかりません")
        print(f"  利用可能なクラス: {sorted(pd.read_excel(path, sheet_name=sheet)['class'].unique())}")
        sys.exit(1)

    # FamilyKey レベルに集約（ラベルあり・なしをまとめる）
    family = df.drop_duplicates(
        subset=["cellline", "treatment", "class", "FamilyKey"]
    ).copy()
    family = family.rename(columns={
        "FamilyKey":                 "species",
        "Family_percent_in_class":   "pct",
    })
    family["group"] = family["cellline"].apply(
        lambda c: "Sensitive" if c in SENSITIVE else "Resistant"
    )

    n_species  = family["species"].nunique()
    n_cells    = family["cellline"].nunique()
    n_trts     = family["treatment"].nunique()
    print(f"[{lipid_class}] 読み込み完了: {n_species} families × "
          f"{n_cells} cell lines × {n_trts} treatments")
    return family


def out_dir(lipid_class, experiment="isotope"):
    d = os.path.join(BASE_RESULTS, experiment, lipid_class)
    os.makedirs(d, exist_ok=True)
    return d


# ---- Figure 1: Stacked bar (Resistant vs Sensitive mean) ----
def plot_stacked_bar(df, lipid_class, treatment, outdir):
    sub = df[df["treatment"] == treatment]
    if len(sub) == 0:
        return

    summary = sub.groupby(["group", "species"])["pct"].mean().reset_index()
    pivot = summary.pivot(index="species", columns="group", values="pct").fillna(0)

    for g in ["Resistant", "Sensitive"]:
        if g not in pivot.columns:
            pivot[g] = 0.0
    pivot = pivot[["Resistant", "Sensitive"]]
    pivot = pivot.sort_values("Sensitive", ascending=False)

    max_pct = pivot.max(axis=1)
    major   = pivot[max_pct >= OTHERS_THRESHOLD]
    minor   = pivot[max_pct <  OTHERS_THRESHOLD]
    if len(minor) > 0:
        pivot_plot = pd.concat([major, minor.sum(axis=0).rename("Others").to_frame().T])
    else:
        pivot_plot = major

    n = len(pivot_plot)
    colors = list(cm.tab20(np.linspace(0, 1, n - (1 if len(minor) > 0 else 0))))
    if len(minor) > 0:
        colors.append((0.75, 0.75, 0.75, 1.0))

    groups_order = ["Resistant", "Sensitive"]
    fig, ax = plt.subplots(figsize=(5, 6))
    bottoms = np.zeros(2)
    for i, sp in enumerate(pivot_plot.index):
        vals = [pivot_plot.loc[sp, g] for g in groups_order]
        ax.bar(groups_order, vals, bottom=bottoms, color=colors[i], label=sp, width=0.5)
        bottoms += np.array(vals)

    ax.set_ylabel(f"mol% of total {lipid_class}", fontsize=12)
    ax.set_title(f"{lipid_class} Composition\n(Resistant n=2 vs Sensitive n=6, {treatment})",
                 fontsize=11)
    ax.set_ylim(0, 105)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8,
              frameon=False, title=f"Species (≥{OTHERS_THRESHOLD}%)", title_fontsize=9)
    plt.tight_layout()

    tag  = TREATMENT_TAGS.get(treatment, treatment.replace(":", "").replace("+", "_"))
    stem = os.path.join(outdir, f"{lipid_class}_isotope_composition_{tag}")
    plt.savefig(stem + ".pdf", bbox_inches="tight")
    plt.savefig(stem + ".png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {stem}.png")


# ---- Figure 2: Dot plot by cell line (top N species) ----
def plot_dot_by_cellline(df, lipid_class, outdir, top_n=10):
    # BSA条件でのSensitive平均が高い上位N種を選出
    bsa = df[df["treatment"] == "BSA"]
    top_species = (
        bsa.groupby("species")["pct"].mean()
        .sort_values(ascending=False)
        .head(top_n)
        .index.tolist()
    )

    avail_trts = [t for t in TREATMENTS if t in df["treatment"].unique()]
    sub = df[df["species"].isin(top_species) & df["treatment"].isin(avail_trts)]

    cellline_order = [c for c in SENSITIVE + RESISTANT if c in df["cellline"].unique()]
    n_sp  = len(top_species)
    n_trt = len(avail_trts)

    fig, axes = plt.subplots(n_sp, n_trt,
                             figsize=(3.0 * n_trt, 2.2 * n_sp),
                             squeeze=False)

    for row_i, sp in enumerate(top_species):
        for col_i, trt in enumerate(avail_trts):
            ax = axes[row_i][col_i]
            sub2 = sub[(sub["species"] == sp) & (sub["treatment"] == trt)]

            for cl_i, cl in enumerate(cellline_order):
                row = sub2[sub2["cellline"] == cl]
                if len(row) == 0:
                    continue
                y     = row["pct"].values[0]
                color = GROUP_COLORS["Resistant"] if cl in RESISTANT else GROUP_COLORS["Sensitive"]
                ax.scatter([cl_i], [y], color=color, s=55, zorder=3, edgecolors="none")

            # グループ平均線
            for grp, cls_in_grp in GROUPS.items():
                idxs = [i for i, c in enumerate(cellline_order) if c in cls_in_grp]
                vals = sub2[sub2["cellline"].isin(cls_in_grp)]["pct"].values
                if len(vals) > 0 and len(idxs) > 0:
                    ax.hlines(vals.mean(), min(idxs) - 0.3, max(idxs) + 0.3,
                              color=GROUP_COLORS[grp], linewidth=1.5, alpha=0.7)

            ax.set_xlim(-0.5, len(cellline_order) - 0.5)
            ax.set_xticks(range(len(cellline_order)))
            ax.set_xticklabels(
                cellline_order, rotation=45, ha="right", fontsize=7
            )
            ax.set_ylim(bottom=0)
            ax.tick_params(axis="y", labelsize=7)

            if col_i == 0:
                ax.set_ylabel(f"{sp}\n(mol%)", fontsize=7)
            if row_i == 0:
                ax.set_title(trt, fontsize=9)

    patches = [
        mpatches.Patch(color=GROUP_COLORS["Sensitive"], label="Sensitive (n=6)"),
        mpatches.Patch(color=GROUP_COLORS["Resistant"], label="Resistant (n=2)"),
    ]
    fig.legend(handles=patches, loc="upper right", fontsize=8, frameon=False)
    fig.suptitle(f"{lipid_class}: Top {top_n} species by cell line",
                 fontsize=12, y=1.005)
    plt.tight_layout()

    stem = os.path.join(outdir, f"{lipid_class}_isotope_dotplot")
    plt.savefig(stem + ".pdf", bbox_inches="tight")
    plt.savefig(stem + ".png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {stem}.png")


# ---- Figure 3: 18:0-d35 取り込み率 (dot plot) ----
def plot_d35_incorporation(lipid_class, experiment, outdir,
                           sheet=DEFAULT_SHEET, top_n=10):
    path = os.path.join(BASE_DATA, experiment, EXCEL_FILE)
    df_raw = pd.read_excel(path, sheet_name=sheet)
    df_raw = df_raw[df_raw["class"] == lipid_class].copy()

    family = df_raw.drop_duplicates(
        subset=["cellline", "treatment", "class", "FamilyKey"]
    ).copy()
    family["group"] = family["cellline"].apply(
        lambda c: "Sensitive" if c in SENSITIVE else "Resistant"
    )

    d35_trts = [t for t in ["18:0-d35", "18:0-d35+18:1-13C5"]
                if t in family["treatment"].unique()]
    if not d35_trts:
        print(f"  [{lipid_class}] d35 treatment なし — スキップ")
        return

    # d35 取り込みが多い上位 N ファミリー
    top_sp = (
        family[family["treatment"].isin(d35_trts)]
        .groupby("FamilyKey")["d35_family_percent_in_class"]
        .mean()
        .sort_values(ascending=False)
        .head(top_n)
        .index.tolist()
    )
    if not top_sp:
        print(f"  [{lipid_class}] d35 incorporation データなし — スキップ")
        return

    sub = family[family["FamilyKey"].isin(top_sp) & family["treatment"].isin(d35_trts)]
    cellline_order = [c for c in SENSITIVE + RESISTANT if c in family["cellline"].unique()]
    n_sp  = len(top_sp)
    n_trt = len(d35_trts)

    fig, axes = plt.subplots(n_sp, n_trt,
                             figsize=(3.0 * n_trt, 2.2 * n_sp),
                             squeeze=False)

    for row_i, sp in enumerate(top_sp):
        for col_i, trt in enumerate(d35_trts):
            ax = axes[row_i][col_i]
            sub2 = sub[(sub["FamilyKey"] == sp) & (sub["treatment"] == trt)]

            for cl_i, cl in enumerate(cellline_order):
                row = sub2[sub2["cellline"] == cl]
                if len(row) == 0:
                    continue
                y     = row["d35_family_percent_in_class"].values[0]
                color = GROUP_COLORS["Resistant"] if cl in RESISTANT else GROUP_COLORS["Sensitive"]
                ax.scatter([cl_i], [y], color=color, s=55, zorder=3, edgecolors="none")

            # グループ平均線
            for grp, cls_in_grp in GROUPS.items():
                idxs = [i for i, c in enumerate(cellline_order) if c in cls_in_grp]
                vals = sub2[sub2["cellline"].isin(cls_in_grp)]["d35_family_percent_in_class"].values
                if len(vals) > 0 and len(idxs) > 0:
                    ax.hlines(vals.mean(), min(idxs) - 0.3, max(idxs) + 0.3,
                              color=GROUP_COLORS[grp], linewidth=1.5, alpha=0.7)

            ax.set_xlim(-0.5, len(cellline_order) - 0.5)
            ax.set_xticks(range(len(cellline_order)))
            ax.set_xticklabels(
                cellline_order, rotation=45, ha="right", fontsize=7
            )
            ax.set_ylim(bottom=0)
            ax.tick_params(axis="y", labelsize=7)

            if col_i == 0:
                ax.set_ylabel(f"{sp}\n(d35 mol%)", fontsize=7)
            if row_i == 0:
                ax.set_title(trt, fontsize=9)

    patches = [
        mpatches.Patch(color=GROUP_COLORS["Sensitive"], label="Sensitive (n=6)"),
        mpatches.Patch(color=GROUP_COLORS["Resistant"], label="Resistant (n=2)"),
    ]
    fig.legend(handles=patches, loc="upper right", fontsize=8, frameon=False)
    fig.suptitle(f"{lipid_class}: 18:0-d35 incorporation (% of family)",
                 fontsize=12, y=1.005)
    plt.tight_layout()

    stem = os.path.join(outdir, f"{lipid_class}_isotope_d35_incorporation")
    plt.savefig(stem + ".pdf", bbox_inches="tight")
    plt.savefig(stem + ".png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {stem}.png")


# ---- Figure 4: Focused species analysis ----
def plot_species_focus(lipid_class, species_name, experiment, outdir,
                       sheet=DEFAULT_SHEET):
    """
    特定の species に着目した focused figure。
    Left : 全治療条件の mol% (lines per cell line, colored by group)
    Right: 18:0-d35 由来の取り込み量（d35 mol%）
    """
    path = os.path.join(BASE_DATA, experiment, EXCEL_FILE)
    df_raw = pd.read_excel(path, sheet_name=sheet)
    df_raw = df_raw[df_raw["class"] == lipid_class].copy()

    family = df_raw.drop_duplicates(
        subset=["cellline", "treatment", "class", "FamilyKey"]
    ).copy()
    family["group"] = family["cellline"].apply(
        lambda c: "Sensitive" if c in SENSITIVE else "Resistant"
    )

    sub = family[family["FamilyKey"] == species_name].copy()
    if len(sub) == 0:
        print(f"  [focus] '{species_name}' が見つかりません")
        return

    avail_trts    = [t for t in TREATMENTS if t in sub["treatment"].unique()]
    d35_trts      = [t for t in ["18:0-d35", "18:0-d35+18:1-13C5"] if t in avail_trts]
    cellline_order = [c for c in SENSITIVE + RESISTANT if c in sub["cellline"].unique()]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # ---- Panel A: mol% across all treatments ----
    ax = axes[0]
    trt_x = {t: i for i, t in enumerate(avail_trts)}

    for cl in cellline_order:
        cl_sub = sub[sub["cellline"] == cl]
        xs = [trt_x[t] for t in avail_trts if t in cl_sub["treatment"].values]
        ys = [float(cl_sub.loc[cl_sub["treatment"] == t, "Family_percent_in_class"].iloc[0])
              for t in avail_trts if t in cl_sub["treatment"].values]
        color = GROUP_COLORS["Resistant"] if cl in RESISTANT else GROUP_COLORS["Sensitive"]
        ax.plot(xs, ys, "o-", color=color, alpha=0.75,
                markersize=8, linewidth=1.3)
        # 最後の点にcell line名をアノテーション
        ax.annotate(cl, xy=(xs[-1], ys[-1]),
                    xytext=(5, 0), textcoords="offset points",
                    fontsize=7.5, color=color, va="center")

    ax.set_xticks(range(len(avail_trts)))
    ax.set_xticklabels(avail_trts, rotation=25, ha="right", fontsize=9)
    ax.set_ylabel(f"mol% of total {lipid_class}", fontsize=11)
    ax.set_title(f"{species_name}\nComposition (all treatments)", fontsize=11)
    ax.set_ylim(bottom=0)
    ax.set_xlim(-0.3, len(avail_trts) - 0.3)

    # ---- Panel B: d35 取り込み量 ----
    ax = axes[1]
    if d35_trts:
        d35_x = {t: i for i, t in enumerate(d35_trts)}

        for cl in cellline_order:
            cl_sub = sub[sub["cellline"] == cl]
            xs = [d35_x[t] for t in d35_trts if t in cl_sub["treatment"].values]
            ys = [float(cl_sub.loc[cl_sub["treatment"] == t,
                                    "d35_family_percent_in_class"].iloc[0])
                  for t in d35_trts if t in cl_sub["treatment"].values]
            color = GROUP_COLORS["Resistant"] if cl in RESISTANT else GROUP_COLORS["Sensitive"]
            ax.plot(xs, ys, "o-", color=color, alpha=0.75,
                    markersize=8, linewidth=1.3)
            ax.annotate(cl, xy=(xs[-1], ys[-1]),
                        xytext=(5, 0), textcoords="offset points",
                        fontsize=7.5, color=color, va="center")

        ax.set_xticks(range(len(d35_trts)))
        ax.set_xticklabels(d35_trts, rotation=25, ha="right", fontsize=9)
        ax.set_ylabel(f"18:0-d35 derived\n(mol% of total {lipid_class})", fontsize=11)
        ax.set_title(f"{species_name}\n18:0-d35 incorporation", fontsize=11)
        ax.set_ylim(bottom=0)
        ax.set_xlim(-0.3, len(d35_trts) - 0.3)
    else:
        ax.set_visible(False)

    # 凡例
    handles = [
        plt.Line2D([0], [0], marker="o", color=GROUP_COLORS["Sensitive"],
                   label="Sensitive (n=6)", linewidth=1.3, markersize=8),
        plt.Line2D([0], [0], marker="o", color=GROUP_COLORS["Resistant"],
                   label="Resistant (n=2)", linewidth=1.3, markersize=8),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2,
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, 1.04))

    fig.suptitle(
        f"{lipid_class}: {species_name}  —  isotope experiment",
        fontsize=12, y=1.08
    )
    plt.tight_layout()

    sp_tag = (species_name.replace("(", "").replace(")", "")
              .replace(":", "").replace("_", "-"))
    stem = os.path.join(outdir, f"{lipid_class}_isotope_focus_{sp_tag}")
    plt.savefig(stem + ".pdf", bbox_inches="tight")
    plt.savefig(stem + ".png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {stem}.png")


# ---------- メイン ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Isotope lipidomics pipeline")
    parser.add_argument("lipid_class",
                        help="脂質クラス名（PI, PC, PE, PS 等）")
    parser.add_argument("--experiment", default="isotope",
                        help="実験名サブディレクトリ（デフォルト: isotope）")
    parser.add_argument("--sheet", default=DEFAULT_SHEET,
                        help=f"シート名（デフォルト: {DEFAULT_SHEET}）")
    parser.add_argument("--focus",
                        help="注目するspecies名（例: 'PI(18:0_18:0)'）。"
                             "指定した場合はfocused figureのみ生成。")
    args = parser.parse_args()

    lipid_class = args.lipid_class
    experiment  = args.experiment

    outdir = out_dir(lipid_class, experiment)

    if args.focus:
        # --focus 指定時: focused figure のみ生成
        print(f"\n===== {lipid_class} focused: {args.focus} [{experiment}] =====")
        plot_species_focus(lipid_class, args.focus, experiment, outdir, args.sheet)
        print(f"\n===== 完了: 出力先 {outdir}/ =====\n")
    else:
        # 通常モード: 全図を生成
        print(f"\n===== {lipid_class} 解析開始 [{experiment}] =====")

        df = load_data(lipid_class, experiment, args.sheet)
        avail_trts = df["treatment"].unique()

        print("[1/3] Stacked bar charts (Resistant vs Sensitive) ...")
        for trt in TREATMENTS:
            if trt in avail_trts:
                plot_stacked_bar(df, lipid_class, trt, outdir)

        print("[2/3] Dot plot by cell line ...")
        plot_dot_by_cellline(df, lipid_class, outdir)

        print("[3/3] 18:0-d35 incorporation ...")
        plot_d35_incorporation(lipid_class, experiment, outdir, args.sheet)

        print(f"\n===== 完了: 出力先 {outdir}/ =====\n")
