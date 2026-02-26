# Lipidomics Composition Analysis

## ディレクトリ構造

```
claudecode_test/
├── README.md
├── pipeline.py          # 解析パイプライン本体
├── stopwatch.html
├── data/                # 入力データ（CSVのみ管理、xlsxはgitignore）
│   ├── SM_pct_individual.csv
│   └── SM_pct_summary.csv
└── results/             # 解析出力（クラスごとにサブフォルダ）
    └── SM/
        ├── *_composition_*.png   # Stacked bar charts
        ├── *_volcano_*.png       # Volcano plots
        └── *_stats_*.csv         # 統計検定結果
```

## 使い方

### 既存クラスを再解析する
```bash
python pipeline.py SM
```

### 新しいクラスを解析する
1. `data/` に `{クラス名}_pct_individual.csv` を置く
2. 以下を実行する：
```bash
python pipeline.py PC
python pipeline.py PE
python pipeline.py PI
python pipeline.py PS
python pipeline.py Cer
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
