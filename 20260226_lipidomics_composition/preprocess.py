"""
Lipidomics Preprocessing: Excel → per-class CSV

Usage:
    # 全クラスを一括処理
    python preprocess.py data/Summary_edit_20250226.xlsx

    # 特定クラスのみ
    python preprocess.py data/Summary_edit_20250226.xlsx SM Cer PC

    # シート名を指定（デフォルト: all_3_重複除去）
    python preprocess.py data/NewExperiment.xlsx --sheet all_3_重複除去

Required files:
    data/sample_metadata.csv  ← サンプル名とメタ情報の対応表
                                 (sample, cellline, sensitivity, group, treatment)

Output:
    data/{CLASS}_pct_individual.csv  ← pipeline.py の入力ファイル
"""

import sys
import os
import argparse
import pandas as pd

DATA_DIR      = "data"
META_FILE     = os.path.join(DATA_DIR, "sample_metadata.csv")
DEFAULT_SHEET = "all_3_重複除去"


def load_metadata():
    if not os.path.exists(META_FILE):
        print(f"[ERROR] メタデータファイルが見つかりません: {META_FILE}")
        print("  sample_metadata.csv を data/ フォルダに置いてください。")
        print("  必要な列: sample, cellline, sensitivity, group, treatment")
        sys.exit(1)
    return pd.read_csv(META_FILE)


def process_class(df, lipid_class, meta, outdir):
    sub = df[df["class"] == lipid_class].copy()
    if len(sub) == 0:
        print(f"  [{lipid_class}] データなし — スキップ")
        return

    # Area列を取得
    area_cols = [c for c in df.columns if str(c).startswith("Area[")]

    # long形式に変換
    long = sub[["species"] + area_cols].melt(
        id_vars="species", var_name="sample", value_name="area"
    )

    # mol% を計算
    totals = long.groupby("sample")["area"].sum().rename("total")
    long = long.join(totals, on="sample")
    long["pct"] = long["area"] / long["total"] * 100

    # メタ情報を結合
    result = long.merge(meta, on="sample", how="inner")
    result = result[["species", "sample", "cellline",
                     "sensitivity", "group", "treatment", "pct"]]

    n_samples  = result["sample"].nunique()
    n_species  = result["species"].nunique()
    out_path   = os.path.join(outdir, f"{lipid_class}_pct_individual.csv")
    result.to_csv(out_path, index=False)
    print(f"  [{lipid_class}] {n_species} species × {n_samples} samples → {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Lipidomics preprocessing")
    parser.add_argument("excel",  help="Excelファイルのパス")
    parser.add_argument("classes", nargs="*",
                        help="処理するクラス名（省略時は全クラス）")
    parser.add_argument("--sheet", default=DEFAULT_SHEET,
                        help=f"シート名（デフォルト: {DEFAULT_SHEET}）")
    args = parser.parse_args()

    # Excel読み込み
    print(f"\nExcel読み込み中: {args.excel}  シート: {args.sheet}")
    if not os.path.exists(args.excel):
        print(f"[ERROR] ファイルが見つかりません: {args.excel}")
        sys.exit(1)

    df = pd.read_excel(args.excel, sheet_name=args.sheet)
    df = df[df["class"].notna()]   # class列がNaNの行を除外
    print(f"  {len(df)} 行 読み込み完了")

    # 対象クラスを決定
    all_classes = sorted(df["class"].unique())
    targets = args.classes if args.classes else all_classes
    print(f"  処理対象クラス: {targets}")

    # メタデータ読み込み
    meta = load_metadata()
    print(f"  メタデータ: {len(meta)} サンプル")

    # クラスごとにCSV出力
    print("\n--- CSV生成 ---")
    for cls in targets:
        process_class(df, cls, meta, DATA_DIR)

    print("\n完了。次は pipeline.py を実行してください:")
    for cls in targets:
        print(f"  python pipeline.py {cls}")


if __name__ == "__main__":
    main()
