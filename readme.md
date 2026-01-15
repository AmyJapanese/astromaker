# astromaker

Pythonで動作する西洋占星術ホロスコープ作成ツールです。  
出生日時と出生地から、惑星・ハウス・軸・各種計算点を算出し、  
数値表示およびチャート描画を行います。

## 特徴
- Swiss Ephemeris（pyswisseph）による高精度計算
- 惑星・小惑星・感受点を同一ロジックで扱う設計
- ハウス・ASC / DSC / MC / IC / Vertex 対応
- Python + matplotlib によるチャート描画
- ephem 依存なし

---

## 動作環境
- Python 3.10 以降（3.11 推奨）
- Windows 環境を想定（他OSは未検証）

---

## インストール

### 必要なライブラリ
以下のライブラリが必要です。

- `swisseph`（pyswisseph）
- `matplotlib`
- `zoneinfo`（Python標準）
- `geopy`
- `timezonefinder`

```bash
pip install pyswisseph matplotlib geopy timezonefinder
````

### pyswisseph のインストールについて

環境によっては `pip install pyswisseph` 実行時に
C++ コンパイラが必要になる場合があります。

その場合は **Microsoft Visual Studio C++ Build Tools** をインストールしてください。

* MSVC v143
* Windows 10 / 11 SDK

を含めれば十分です。

---

## 使い方

```bash
python astromaker.py
```

起動後、出生日時・出生地を入力するとホロスコープを生成します。

### 小惑星について

小惑星を正確に出すには ephemeris ファイル（seas_18.se1 など）を SWE_EPH_PATH に置くと安定します。

ただしSwiss ephemerisの仕様上、絶対パスに日本語を含めないでください。

---

## 注意

本ツールでは、グランドクロスは構造上Tスクエア4つとしても検出されます。

これは誤検出ではなく、グランドクロスが複数のTスクエアを内包するという力学をそのまま出力しているためです。

表示が冗長になる場合がありますが、仕様です。ご了承ください（そして我慢してください）。

---

## ライセンス

MIT License