# Lipidomics Composition Analysis — Results Summary

## 解析概要

- **目的**: C18:0処理に対するリピドミクス応答のクラス別コンポジション解析
- **比較軸**:
  - A) Resistant vs Sensitive（感受性の違い）
  - B) Ctr vs C18:0（処理効果）
- **統計検定**: Mann-Whitney U test、FDR補正: Benjamini-Hochberg法（q < 0.2 を有意とした）
- **サンプル数**: Resistant 7〜8株、Sensitive 13〜19株
- **解析クラス**: SM, Cer, PC, PE, PI, PS

---

## クラス別結果

### SM（Sphingomyelin）| 83 species

**A) Resistant vs Sensitive**: FDR有意差なし（Ctr・C18:0ともに）

**B) C18:0処理効果**:

| グループ | FDR q<0.2 | 傾向 |
|---------|-----------|------|
| Resistant | **0種** | 変化なし |
| Sensitive | **16種**（↑4, ↓12） | 長鎖飽和SMが増加 |

- Sensitive で増加した主要種: `SM(d38:1)` (log2FC=+1.97), `SM(d36:0)` (+1.35), `SM(d38:0)` (+1.33), `SM(d18:1_18:0)` (+0.98)
- **→ C18:0処理への応答はSensitiveに選択的**

---

### Cer（Ceramide）| 76 species

**A) Resistant vs Sensitive**: FDR有意差なし（Ctr・C18:0ともに）

**B) C18:0処理効果**:

| グループ | FDR q<0.2 | 傾向 |
|---------|-----------|------|
| Resistant | **3種**（↑0, ↓3） | 一部マイナー種が低下 |
| Sensitive | **19種**（↑15, ↓4） | 飽和・長鎖Cerが増加 |

- Sensitive で増加した主要種: `Cer(t18:1_14:0)` (log2FC=+3.46), `Cer(m20:2_18:0)` (+2.64), `Cer(d20:0_18:0)` (+1.91), `Cer(d18:0_18:0)` (+1.72)
- **→ SMと同様、飽和長鎖種の増加がSensitiveで顕著**

---

### PC（Phosphatidylcholine）| 476 species

**A) Resistant vs Sensitive**: FDR有意差なし（Ctr・C18:0ともに）

**B) C18:0処理効果**:

| グループ | FDR q<0.2 | 傾向 |
|---------|-----------|------|
| Resistant | **72種**（↑21, ↓51） | 多種が低下 |
| Sensitive | **222種**（↑52, ↓170） | より広範囲に低下 |

- **→ 両グループともに応答するが、Sensitiveの方が規模が大きい**
- PCは全体的に低下傾向（不飽和PC減少）

---

### PE（Phosphatidylethanolamine）| 241 species

**A) Resistant vs Sensitive**:
- Ctr: FDR有意差なし
- C18:0: **27種**が有意（主にResistantで低値）

**B) C18:0処理効果**:

| グループ | FDR q<0.2 | 傾向 |
|---------|-----------|------|
| Resistant | **9種**（↑6, ↓3） | 限定的な応答 |
| Sensitive | **117種**（↑27, ↓90） | 大規模な変化 |

- Sensitive で増加した主要種: `PE(23:0_16:1)` (log2FC=+1.80), `PE(18:1_18:0)` (+0.54)
- **→ SM・Cerと同様のSensitive選択的応答パターン**

---

### PI（Phosphatidylinositol）| 78 species

**A) Resistant vs Sensitive**: FDR有意差なし（Ctr）、1種（C18:0）

**B) C18:0処理効果**:

| グループ | FDR q<0.2 | 傾向 |
|---------|-----------|------|
| Resistant | **34種**（↑3, ↓31） | PI全体が低下 |
| Sensitive | **50種**（↑5, ↓45） | PI全体が低下 |

- **→ 両グループともに同方向（低下）で応答。Resistant/Sensitiveの差が小さい**
- PIはC18:0処理で全般的に減少するクラス

---

### PS（Phosphatidylserine）| 55 species

**A) Resistant vs Sensitive**:
- Ctr: FDR有意差なし
- C18:0: **6種**が有意（q < 0.2）

**B) C18:0処理効果**:

| グループ | FDR q<0.2 | 傾向 |
|---------|-----------|------|
| Resistant | **28種**（↑2, ↓26） | 大半が低下 |
| Sensitive | **34種**（↑4, ↓30） | 大半が低下 |

- 両グループで `PS(18:0_18:1)` が顕著に増加（共通応答）
- **→ PIと同様に両グループが応答。PS(18:0_18:1)の上昇は共通の特徴**

---

## 総括

### C18:0処理への応答パターン

| クラス | Resistant応答 | Sensitive応答 | パターン |
|--------|:---:|:---:|---------|
| SM | ✗ | ✓ | **Sensitive選択的** |
| Cer | △（3種） | ✓ | **Sensitive選択的** |
| PE | △（9種） | ✓ | **Sensitive選択的** |
| PC | ✓ | ✓✓ | 両者応答・Sensitive優位 |
| PI | ✓ | ✓ | 両者応答・方向性一致（低下） |
| PS | ✓ | ✓ | 両者応答・方向性一致（低下） |

### 主要な知見

1. **SM・Cer・PEで一貫した「Sensitive選択的応答」パターン**
   - C18:0処理後に飽和・長鎖種が増加するのはSensitiveのみ
   - Resistantではこの代謝リモデリングが起きない

2. **PI・PSはResistant/Sensitiveを問わず広範な低下**
   - C18:0処理の非選択的な効果の可能性

3. **Resistant vs Sensitive の基礎的なコンポジションの差は小さい**
   - 無処理時点ではほとんどのクラスで有意差なし
   - 差が現れるのはC18:0処理後（PE・PS）
