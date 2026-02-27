# Lipidomics Composition Analysis

## ディレクトリ構造

```
20260226_lipidomics_composition/
├── README.md
├── preprocess.py        # Step 1: Excel → per-class CSV
├── pipeline.py          # Step 2: コンポジション解析・図の生成
├── stopwatch.html
├── data/
│   ├── sample_metadata.csv        # サンプル情報（実験ごとに差し替え）
│   ├── Summary_edit_20250226.xlsx # 元データ（gitignore）
│   ├── SM_pct_individual.csv      # pipeline.py の入力
│   ├── Cer_pct_individual.csv
│   ├── PC_pct_individual.csv
│   └── ...
└── results/                       # 解析出力（クラスごとにサブフォルダ）
    ├── SM/
    ├── Cer/
    └── ...
```

## ワークフロー

### Step 1: ExcelからCSVを生成（実験ごとに1回）

```bash
# 全クラスを一括処理
python preprocess.py data/Summary_edit_20250226.xlsx

# 特定クラスのみ
python preprocess.py data/Summary_edit_20250226.xlsx SM Cer PC PE PI PS

# シート名が異なる場合
python preprocess.py data/NewExperiment.xlsx --sheet シート名
```

**前提ファイル：** `data/sample_metadata.csv`

| 列名 | 内容 | 例 |
|------|------|----|
| sample | サンプル列名 | Area[s11-1] |
| cellline | 細胞株名 | A549 |
| sensitivity | 感受性グループ | Res / C18 / C18-C16 |
| group | Resistant / Sensitive | Resistant |
| treatment | 処理条件 | Ctr / C18:0 / C18:0+C18:1 |

> 新しい実験では `sample_metadata.csv` を差し替えるだけでOK

### Step 2: コンポジション解析の実行

```bash
python pipeline.py SM
python pipeline.py Cer
python pipeline.py PC   # PC, PE, PI, PS も同様
```

→ `results/{クラス名}/` に figure (PNG) と stats (CSV) が自動生成される

## 解析内容

| Figure | 内容 |
|--------|------|
| `*_composition_Ctr.png` | Stacked bar — Resistant vs Sensitive（無処理） |
| `*_composition_C180.png` | Stacked bar — Resistant vs Sensitive（C18:0処理） |
| `*_volcano_ResVsSen.png` | Volcano plot — Resistant vs Sensitive（Ctr & C18:0） |
| `*_volcano_treatment.png` | Volcano plot — Ctr vs C18:0（Resistant & Sensitive） |

統計検定: Mann-Whitney U test、FDR補正: Benjamini-Hochberg法

## Git / GitHub

- リポジトリ: `git@github.com:bathjack/bathjack.git`（SSH接続）
- xlsxファイル・PDFは `.gitignore` で除外済み
