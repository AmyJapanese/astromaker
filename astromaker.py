import tkinter as tk
from tkinter import ttk

import math
import os
from itertools import combinations

import swisseph as swe

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from dataclasses import dataclass
from datetime import datetime, timezone
from zoneinfo import ZoneInfo

from geopy.geocoders import Nominatim
from timezonefinder import TimezoneFinder

# ============================================================
# Swiss Ephemeris: ephe path（小惑星ファイルを置く場所）
#  - 置いてなくても主要10天体は計算できることが多い（内蔵/自動フォールバック）
#  - 小惑星（Ceres等）は ephe ファイルが無いと失敗することがある
# ============================================================
def init_swe_ephe_path() -> str:
    """
    どこに ephemeris ファイルを置くか：
      1) 環境変数 SWEPH_PATH
      2) スクリプトと同じ場所の ./ephe
      3) スクリプトと同じ場所（最後の手段）
    """
    base = os.path.dirname(os.path.abspath(__file__))
    ephe_path = os.environ.get("SWEPH_PATH")
    if ephe_path and os.path.isdir(ephe_path):
        swe.set_ephe_path(ephe_path)
        return ephe_path

    ephe_dir = os.path.join(base, "ephe")
    if os.path.isdir(ephe_dir):
        swe.set_ephe_path(ephe_dir)
        return ephe_dir

    # 何も無ければ base を渡しておく（ファイルがあれば拾える）
    swe.set_ephe_path(base)
    return base

SWE_EPH_PATH = init_swe_ephe_path()

# ============================================================
# 小物：角度の正規化
# ============================================================
def norm360(deg: float) -> float:
    return deg % 360.0

# ============================================================
# Julian Day（UTC）
# ============================================================
def julian_day_ut(dt_utc: datetime) -> float:
    """
    Swiss Ephemeris と相性の良い作り：swe.julday を使う
    dt_utc: timezone-aware (UTC) を想定
    """
    if dt_utc.tzinfo is None:
        raise ValueError("dt_utc must be timezone-aware (UTC).")
    dt_utc = dt_utc.astimezone(timezone.utc)
    hour = dt_utc.hour + dt_utc.minute / 60.0 + (dt_utc.second + dt_utc.microsecond / 1e6) / 3600.0
    return swe.julday(dt_utc.year, dt_utc.month, dt_utc.day, hour, swe.GREG_CAL)

# ============================================================
# 星座
# ============================================================
SIGN_NAMES = ["Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
              "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces"]
SIGN_SYMBOLS = ["♈︎","♉︎","♊︎","♋︎","♌︎","♍︎","♎︎","♏︎","♐︎","♑︎","♒︎","♓︎"]

def deg_to_sign_index(lon_deg: float) -> int:
    return int((lon_deg % 360) // 30)

def deg_to_sign_pos(lon_deg: float):
    lon = lon_deg % 360
    sign_idx = deg_to_sign_index(lon)
    in_sign = lon - sign_idx * 30
    deg_int = int(in_sign)
    minute = int(round((in_sign - deg_int) * 60))
    if minute == 60:
        deg_int += 1
        minute = 0
        if deg_int == 30:
            deg_int = 0
            sign_idx = (sign_idx + 1) % 12
    return sign_idx, deg_int, minute

def deg_to_sign_string(lon_deg: float) -> str:
    s, d, m = deg_to_sign_pos(lon_deg)
    return f"{SIGN_SYMBOLS[s]} {d:02d}°{m:02d}' ({SIGN_NAMES[s]})"

def lon_to_house_equal(lon_deg: float, ac_deg: float) -> int:
    rel = (lon_deg - ac_deg) % 360
    return int(rel // 30) + 1

def choose_house_system(lat: float) -> tuple[str, str]:
    """
    緯度に応じてハウスシステムを選択
    戻り値: (hsys_code, label)
    """
    if abs(lat) >= 60.0:
        return "O", "Porphyry (auto: high latitude)"
    else:
        return "P", "Placidus"

# ============================================================
# アスペクト（メジャー/マイナー）
# ============================================================
MAJOR_ASPECTS = {
    "Conjunction": 0,
    "Sextile": 60,
    "Square": 90,
    "Trine": 120,
    "Opposition": 180,
}

MINOR_ASPECTS = {
    "Semisextile": 30,
    "Semisquare": 45,
    "Sesquiquadrate": 135,
    "Quincunx": 150,
    "Quintile": 72,
    "Biquintile": 144,
}

ASPECT_COLORS = {
    "Conjunction": "gray",
    "Sextile": "dodgerblue",
    "Square": "red",
    "Trine": "limegreen",
    "Opposition": "orange",

    "Semisextile": "lightskyblue",
    "Semisquare": "hotpink",
    "Sesquiquadrate": "violet",
    "Quincunx": "gold",
    "Quintile": "turquoise",
    "Biquintile": "teal",
}

ASPECT_SYMBOLS = {
    "Conjunction": "☌",
    "Sextile": "⚹",
    "Square": "□",
    "Trine": "△",
    "Opposition": "☍",
    "Semisextile": "⚺",
    "Semisquare": "∠",
    "Sesquiquadrate": "⚼",
    "Quincunx": "⚻",
    "Quintile": "Q",
    "Biquintile": "bQ",
}

DEFAULT_MAJOR_ORB_DEG = 6.0
DEFAULT_MINOR_ORB_DEG = 2.0

# ============================================================
# 除外ルール：特定ペアの Opposition（180°）は「無いもの」として扱う
#  - ASC & DSC
#  - MC  & IC
#  - NorthNode & SouthNode
#  - Lilith & Selena
# ============================================================
EXCLUDED_OPPOSITION_PAIRS = {
    frozenset(("ASC", "DSC")),
    frozenset(("MC", "IC")),
    frozenset(("NorthNode", "SouthNode")),
    frozenset(("Lilith", "Selena")),
}

AXES_MAIN = {"ASC", "DSC", "MC", "IC"}

def is_excluded_opposition_pair(n1: str, n2: str) -> bool:
    return frozenset((n1, n2)) in EXCLUDED_OPPOSITION_PAIRS

# ============================================================
# チャートパターン検出（コンソール用）
# ============================================================

def _sep(a: float, b: float) -> float:
    return angular_separation_deg(a, b)

def _hit(a: float, b: float, target: float, orb: float):
    """return (ok, delta, sep)"""
    s = _sep(a, b)
    d = abs(s - target)
    return (d <= orb), d, s

def _get_positions(longitudes: dict, enabled_names: list[str]) -> dict[str, float]:
    """None を落として position dict を作る"""
    out = {}
    for n in enabled_names:
        v = longitudes.get(n)
        if v is not None:
            out[n] = float(v) % 360.0
    return out

def _dedup_tuples(items):
    """list[tuple] をソート＋set で重複排除"""
    seen = set()
    out = []
    for t in items:
        key = tuple(t)
        if key in seen:
            continue
        seen.add(key)
        out.append(t)
    return out

def detect_grand_trines(pos: dict[str, float], orb_trine: float):
    hits = []
    names = list(pos.keys())
    for a, b, c in combinations(names, 3):
        ok1, d1, _ = _hit(pos[a], pos[b], 120, orb_trine)
        ok2, d2, _ = _hit(pos[b], pos[c], 120, orb_trine)
        ok3, d3, _ = _hit(pos[c], pos[a], 120, orb_trine)
        if ok1 and ok2 and ok3:
            score = max(d1, d2, d3)
            hits.append((score, tuple(sorted((a, b, c)))))
    hits.sort(key=lambda x: x[0])
    return [t for _, t in hits]

def detect_t_squares(pos: dict[str, float], orb_opp: float, orb_sq: float):
    """
    A-B opposition, C squares both => (baseA, baseB, apexC)
    """
    hits = []
    names = list(pos.keys())
    for a, b in combinations(names, 2):
        ok_opp, d_opp, _ = _hit(pos[a], pos[b], 180, orb_opp)
        if ok_opp and is_excluded_opposition_pair(a, b):
            ok_opp = False
        if not ok_opp:
            continue
        for c in names:
            if c == a or c == b:
                continue
            ok1, d1, _ = _hit(pos[a], pos[c], 90, orb_sq)
            ok2, d2, _ = _hit(pos[b], pos[c], 90, orb_sq)
            if ok1 and ok2:
                score = max(d_opp, d1, d2)
                base = tuple(sorted((a, b)))
                hits.append((score, base[0], base[1], c))
    hits.sort(key=lambda x: x[0])
    # 重複（base同じ＆apex同じ）を消す
    seen = set()
    out = []
    for score, a, b, c in hits:
        key = (a, b, c)
        if key in seen:
            continue
        seen.add(key)
        out.append((a, b, c))
    return out

def detect_grand_crosses(pos: dict[str, float], orb_opp: float, orb_sq: float):
    """
    4天体で：opposition 2本 + square 4本（合計6ペア）
    """
    hits = []
    names = list(pos.keys())
    for a, b, c, d in combinations(names, 4):
        pair = [(a, b), (a, c), (a, d), (b, c), (b, d), (c, d)]
        opp = 0
        sq = 0
        worst = 0.0
        ok = True
        for x, y in pair:
            ok_opp, d_opp, _ = _hit(pos[x], pos[y], 180, orb_opp)
            if ok_opp and is_excluded_opposition_pair(x, y):
                ok_opp = False
            if ok_opp and is_excluded_opposition_pair(x, y):
                ok_opp = False
            if ok_opp and is_excluded_opposition_pair(x, y):
                ok_opp = False
            ok_sq, d_sq, _ = _hit(pos[x], pos[y], 90, orb_sq)
            if ok_opp:
                opp += 1
                worst = max(worst, d_opp)
            elif ok_sq:
                sq += 1
                worst = max(worst, d_sq)
            else:
                ok = False
                break
        if ok and opp == 2 and sq == 4:
            hits.append((worst, tuple(sorted((a, b, c, d)))))
    hits.sort(key=lambda x: x[0])
    return [t for _, t in hits]

def detect_kites(pos: dict[str, float], grand_trines: list[tuple[str, str, str]],
                 orb_opp: float, orb_sex: float):
    """
    Kite: Grand Trine (A,B,C) + D opposite one vertex (focus) AND D sextile the other two
    => (A,B,C,D, focus)
    """
    hits = []
    names_all = set(pos.keys())
    for a, b, c in grand_trines:
        tri = (a, b, c)
        for focus in tri:
            others = [x for x in tri if x != focus]
            for d in names_all:
                if d in tri:
                    continue
                ok_opp, d_opp, _ = _hit(pos[d], pos[focus], 180, orb_opp)
                if ok_opp and is_excluded_opposition_pair(d, focus):
                    ok_opp = False
                if not ok_opp:
                    continue
                ok_s1, d1, _ = _hit(pos[d], pos[others[0]], 60, orb_sex)
                ok_s2, d2, _ = _hit(pos[d], pos[others[1]], 60, orb_sex)
                if ok_s1 and ok_s2:
                    score = max(d_opp, d1, d2)
                    hits.append((score, tuple(sorted(tri)), d, focus))
    hits.sort(key=lambda x: x[0])

    out = []
    seen = set()
    for score, tri_sorted, d, focus in hits:
        key = (tri_sorted, d, focus)
        if key in seen:
            continue
        seen.add(key)
        out.append((*tri_sorted, d, focus))
    return out

def detect_mystic_rectangles(pos: dict[str, float], orb_opp: float, orb_trine: float, orb_sex: float):
    """
    Mystic Rectangle: 4天体で opposition 2本 + trine 2本 + sextile 2本
    """
    hits = []
    names = list(pos.keys())
    for a, b, c, d in combinations(names, 4):
        pair = [(a, b), (a, c), (a, d), (b, c), (b, d), (c, d)]
        cnt = {"opp": 0, "tri": 0, "sex": 0}
        worst = 0.0
        ok = True
        for x, y in pair:
            ok_opp, d_opp, _ = _hit(pos[x], pos[y], 180, orb_opp)
            ok_tri, d_tri, _ = _hit(pos[x], pos[y], 120, orb_trine)
            ok_sex, d_sex, _ = _hit(pos[x], pos[y], 60, orb_sex)
            if ok_opp:
                cnt["opp"] += 1
                worst = max(worst, d_opp)
            elif ok_tri:
                cnt["tri"] += 1
                worst = max(worst, d_tri)
            elif ok_sex:
                cnt["sex"] += 1
                worst = max(worst, d_sex)
            else:
                ok = False
                break
        if ok and cnt["opp"] == 2 and cnt["tri"] == 2 and cnt["sex"] == 2:
            hits.append((worst, tuple(sorted((a, b, c, d)))))
    hits.sort(key=lambda x: x[0])
    return [t for _, t in hits]

def detect_yods(pos: dict[str, float], orb_sex: float, orb_quinc: float):
    """
    Yod: A-B sextile + C quincunx to both => (A,B, apexC)
    """
    hits = []
    names = list(pos.keys())
    for a, b in combinations(names, 2):
        ok_sex, d_sex, _ = _hit(pos[a], pos[b], 60, orb_sex)
        if not ok_sex:
            continue
        for c in names:
            if c == a or c == b:
                continue
            ok1, d1, _ = _hit(pos[a], pos[c], 150, orb_quinc)
            ok2, d2, _ = _hit(pos[b], pos[c], 150, orb_quinc)
            if ok1 and ok2:
                score = max(d_sex, d1, d2)
                base = tuple(sorted((a, b)))
                hits.append((score, base[0], base[1], c))
    hits.sort(key=lambda x: x[0])

    out = []
    seen = set()
    for score, a, b, c in hits:
        key = (a, b, c)
        if key in seen:
            continue
        seen.add(key)
        out.append((a, b, c))
    return out

def detect_boomerangs(pos: dict[str, float], yods: list[tuple[str, str, str]], orb_opp: float):
    """
    Boomerang: Yod + 4th planet opposite apex
    => (A,B, apexC, D_opposite_apex)
    """
    hits = []
    names = list(pos.keys())
    for a, b, apex in yods:
        for d in names:
            if d in (a, b, apex):
                continue
            ok_opp, d_opp, _ = _hit(pos[d], pos[apex], 180, orb_opp)
            if ok_opp and is_excluded_opposition_pair(d, apex):
                ok_opp = False
            if ok_opp:
                hits.append((d_opp, a, b, apex, d))
    hits.sort(key=lambda x: x[0])

    out = []
    seen = set()
    for score, a, b, apex, d in hits:
        key = (tuple(sorted((a, b))), apex, d)
        if key in seen:
            continue
        seen.add(key)
        base = tuple(sorted((a, b)))
        out.append((base[0], base[1], apex, d))
    return out

def detect_golden_yods(pos: dict[str, float], orb_quint: float, orb_biquint: float):
    """
    Golden Yod / Quintile Yod（よくある定義）:
      A-B biquintile(144) + C quintile(72) to both A and B
    => (A,B, apexC)
    """
    hits = []
    names = list(pos.keys())
    for a, b in combinations(names, 2):
        ok_bq, d_bq, _ = _hit(pos[a], pos[b], 144, orb_biquint)
        if not ok_bq:
            continue
        for c in names:
            if c == a or c == b:
                continue
            ok1, d1, _ = _hit(pos[a], pos[c], 72, orb_quint)
            ok2, d2, _ = _hit(pos[b], pos[c], 72, orb_quint)
            if ok1 and ok2:
                score = max(d_bq, d1, d2)
                base = tuple(sorted((a, b)))
                hits.append((score, base[0], base[1], c))
    hits.sort(key=lambda x: x[0])

    out = []
    seen = set()
    for score, a, b, c in hits:
        key = (a, b, c)
        if key in seen:
            continue
        seen.add(key)
        out.append((a, b, c))
    return out

def detect_thors_hammers(pos: dict[str, float], orb_sq: float, orb_sesqui: float):
    """
    Thor's Hammer（よくある定義）:
      A-B square(90) + C sesquiquadrate(135) to both A and B
    => (A,B, apexC)
    """
    hits = []
    names = list(pos.keys())
    for a, b in combinations(names, 2):
        ok_sq, d_sq, _ = _hit(pos[a], pos[b], 90, orb_sq)
        if not ok_sq:
            continue
        for c in names:
            if c == a or c == b:
                continue
            ok1, d1, _ = _hit(pos[a], pos[c], 135, orb_sesqui)
            ok2, d2, _ = _hit(pos[b], pos[c], 135, orb_sesqui)
            if ok1 and ok2:
                score = max(d_sq, d1, d2)
                base = tuple(sorted((a, b)))
                hits.append((score, base[0], base[1], c))
    hits.sort(key=lambda x: x[0])

    out = []
    seen = set()
    for score, a, b, c in hits:
        key = (a, b, c)
        if key in seen:
            continue
        seen.add(key)
        out.append((a, b, c))
    return out

def detect_grand_sextiles(pos: dict[str, float], orb_sex: float, orb_tri: float, orb_opp: float):
    """
    Grand Sextile（スター・オブ・デイビッド）
    理想形のペア数（6天体=15ペア）:
      sextile 6本 / trine 6本 / opposition 3本
    """
    hits = []
    names = list(pos.keys())
    if len(names) < 6:
        return []
    for sext in combinations(names, 6):
        sext = list(sext)
        cnt = {"sex": 0, "tri": 0, "opp": 0}
        worst = 0.0
        ok = True
        for x, y in combinations(sext, 2):
            ok_sex, d_sex, _ = _hit(pos[x], pos[y], 60, orb_sex)
            ok_tri, d_tri, _ = _hit(pos[x], pos[y], 120, orb_tri)
            ok_opp, d_opp, _ = _hit(pos[x], pos[y], 180, orb_opp)
            if ok_sex:
                cnt["sex"] += 1
                worst = max(worst, d_sex)
            elif ok_tri:
                cnt["tri"] += 1
                worst = max(worst, d_tri)
            elif ok_opp:
                cnt["opp"] += 1
                worst = max(worst, d_opp)
            else:
                ok = False
                break
        if ok and cnt["sex"] == 6 and cnt["tri"] == 6 and cnt["opp"] == 3:
            hits.append((worst, tuple(sorted(sext))))
    hits.sort(key=lambda x: x[0])
    # 6天体検出は重いので、同じセットは重複排除
    return [t for _, t in hits]

def detect_stelliums(pos: dict[str, float], min_bodies: int = 3, max_span_deg: float = 10.0):
    """
    Stellium（合の集中）: 円周上で span <= max_span_deg に min_bodies 以上。
    """
    names = list(pos.keys())
    if len(names) < min_bodies:
        return []

    # 0-360 を 2周分にして窓探索（円周対策）
    pts = sorted((pos[n], n) for n in names)
    pts2 = pts + [(lon + 360.0, n) for lon, n in pts]

    hits = []
    # スライディングウィンドウ
    j = 0
    for i in range(len(pts)):
        while j < len(pts2) and (pts2[j][0] - pts2[i][0]) <= max_span_deg:
            j += 1
        window = pts2[i:j]
        # 同じ天体が2周分で混ざるのを排除
        uniq = []
        seen = set()
        for lon, n in window:
            if n in seen:
                continue
            seen.add(n)
            uniq.append(n)
        if len(uniq) >= min_bodies:
            # 正規化したキーで重複排除
            key = tuple(sorted(uniq))
            hits.append(key)

    hits = _dedup_tuples(sorted(set(hits), key=lambda t: (-len(t), t)))
    return hits

def detect_distribution_pattern(pos: dict[str, float]):
    """
    ざっくり配置判定（厳密定義は流派で違うので“目安”）
    - 最大ギャップ max_gap
    - 全体スパン span = 360 - max_gap
    """
    names = list(pos.keys())
    if len(names) < 3:
        return ("(too few bodies)", {})

    lons = sorted((pos[n], n) for n in names)
    gaps = []
    for i in range(len(lons)):
        lon1, _ = lons[i]
        lon2, _ = lons[(i + 1) % len(lons)]
        gap = (lon2 - lon1) % 360.0
        gaps.append(gap)

    max_gap = max(gaps)
    max_i = gaps.index(max_gap)
    span = 360.0 - max_gap

    # 最大ギャップの後から始めると「最小スパン区間」になる
    start_idx = (max_i + 1) % len(lons)
    ordered = lons[start_idx:] + lons[:start_idx]

    info = {
        "span_deg": round(span, 2),
        "max_gap_deg": round(max_gap, 2),
        "bodies": [n for _, n in ordered],
    }

    # 判定（目安）
    # Bundle: 120°以内
    if span <= 120:
        return ("Bundle (集中)", info)

    # Bowl/Bucket/Locomotive: 180°〜240° くらいが多い
    if 120 < span <= 180:
        # Bowl寄り
        # Bucket: 1天体だけ“取っ手”っぽく離れる（隣接ギャップが大きい）
        # ここでは「最大ギャップ以外にも大きいギャップ(>=30°)が1つある」等で雑判定
        big_gaps = sum(1 for g in gaps if g >= 30)
        if big_gaps >= 2:
            return ("Bucket-ish (ボウル＋取っ手っぽい)", info)
        return ("Bowl (半分に偏る)", info)

    if 180 < span <= 240:
        # Locomotive: 240°以内に収まり、空白が 60°以上
        if max_gap >= 60:
            return ("Locomotive (機関車型)", info)
        return ("Bowl-ish (広め)", info)

    # Seesaw: 大きいギャップが複数（2クラスター）
    big_gaps = [g for g in gaps if g >= 60]
    if len(big_gaps) >= 2:
        return ("Seesaw (シーソー/二極)", info)

    # Splash: 全体に散る＆最大ギャップが小さい
    if span > 240 and max_gap < 60:
        return ("Splash (散開)", info)

    return ("(no strong distribution pattern)", info)

def detect_chart_patterns(longitudes: dict, enabled_names: list[str], major_orb: float, minor_orb: float):
    """
    まとめて検出して dict で返す。
    major_orb: メジャー（0/60/90/120/180）
    minor_orb: マイナー（30/45/72/135/144/150）
    """
    pos = _get_positions(longitudes, enabled_names)

    # 何もない時の安全策
    if len(pos) < 3:
        return {"_pos_count": len(pos), "note": "too few bodies"}

    gtr = detect_grand_trines(pos, orb_trine=major_orb)
    tsq = detect_t_squares(pos, orb_opp=major_orb, orb_sq=major_orb)
    gcx = detect_grand_crosses(pos, orb_opp=major_orb, orb_sq=major_orb)
    # 軸4つ(ASC/DSC/MC/IC)だけで作られるグランドクロスは除外
    gcx = [t for t in gcx if set(t) != AXES_MAIN]
    kit = detect_kites(pos, gtr, orb_opp=major_orb, orb_sex=major_orb)
    mrx = detect_mystic_rectangles(pos, orb_opp=major_orb, orb_trine=major_orb, orb_sex=major_orb)

    yod = detect_yods(pos, orb_sex=major_orb, orb_quinc=minor_orb)
    boo = detect_boomerangs(pos, yod, orb_opp=major_orb)

    gyd = detect_golden_yods(pos, orb_quint=minor_orb, orb_biquint=minor_orb)
    thh = detect_thors_hammers(pos, orb_sq=major_orb, orb_sesqui=minor_orb)

    gsx = detect_grand_sextiles(pos, orb_sex=major_orb, orb_tri=major_orb, orb_opp=major_orb)

    # Stellium: span は少し広めでも良いので（好みで調整）
    stellium = detect_stelliums(pos, min_bodies=3, max_span_deg=max(8.0, major_orb * 2))

    dist_name, dist_info = detect_distribution_pattern(pos)

    return {
        "_pos_count": len(pos),
        "GrandTrine": gtr,
        "Kite": kit,
        "TSquare": tsq,
        "GrandCross": gcx,
        "MysticRectangle": mrx,
        "Yod": yod,
        "Boomerang": boo,
        "GoldenYod": gyd,
        "ThorsHammer": thh,
        "GrandSextile": gsx,
        "Stellium": stellium,
        "Distribution": (dist_name, dist_info),
    }

def print_chart_patterns(patterns: dict):
    print("\n=== Chart Patterns (AUTO DETECT) ===")
    if not patterns or patterns.get("_pos_count", 0) < 3:
        print("  (too few bodies)")
        return

    def _print_list(title, items, fmt=None):
        if not items:
            return
        print(f"\n-- {title} ({len(items)}) --")
        for it in items:
            if fmt:
                print("  " + fmt(it))
            else:
                print("  " + ", ".join(it))

    _print_list("Grand Trine", patterns.get("GrandTrine"))
    _print_list("Kite", patterns.get("Kite"),
                fmt=lambda t: f"GT[{', '.join(t[:3])}] + D={t[3]} opposite focus={t[4]}")
    _print_list("T-Square", patterns.get("TSquare"),
                fmt=lambda t: f"base={t[0]}-{t[1]} (opp), apex={t[2]}")
    _print_list("Grand Cross", patterns.get("GrandCross"))
    _print_list("Mystic Rectangle", patterns.get("MysticRectangle"))
    _print_list("Yod", patterns.get("Yod"),
                fmt=lambda t: f"base={t[0]}-{t[1]} (sextile), apex={t[2]} (quincunx×2)")
    _print_list("Boomerang (Yod + opp apex)", patterns.get("Boomerang"),
                fmt=lambda t: f"base={t[0]}-{t[1]}, apex={t[2]}, opp_apex={t[3]}")
    _print_list("Golden Yod (Quintile Yod)", patterns.get("GoldenYod"),
                fmt=lambda t: f"base={t[0]}-{t[1]} (biquintile), apex={t[2]} (quintile×2)")
    _print_list("Thor's Hammer", patterns.get("ThorsHammer"),
                fmt=lambda t: f"base={t[0]}-{t[1]} (square), apex={t[2]} (sesqui×2)")
    _print_list("Grand Sextile (Star of David)", patterns.get("GrandSextile"))

    # Stellium
    stell = patterns.get("Stellium", [])
    if stell:
        print(f"\n-- Stellium (>=3 in tight span) ({len(stell)}) --")
        for t in stell:
            print("  " + ", ".join(t))

    # Distribution
    dist = patterns.get("Distribution")
    if dist:
        name, info = dist
        print("\n-- Distribution (rough) --")
        print(f"  {name}")
        if isinstance(info, dict) and info:
            print(f"  span={info.get('span_deg')}° max_gap={info.get('max_gap_deg')}°")
            # bodies list は長いので、欲しければコメントアウト外してね
            # print("  bodies:", ", ".join(info.get("bodies", [])))

    # 何も出なかったとき
    any_found = False
    for k, v in patterns.items():
        if k.startswith("_") or k == "Distribution":
            continue
        if v:
            any_found = True
            break
    if not any_found and not stell:
        print("  (none found with current orbs)")

# ============================================================
# 場所解決（地名 → 緯度経度 → タイムゾーン）
# ============================================================
@dataclass
class LocationInfo:
    name: str
    lat: float
    lon: float
    tz: str  # IANA timezone

_geolocator = Nominatim(user_agent="astromaker/2.0 (pyswisseph)")
_tzf = TimezoneFinder()
_location_cache: dict[str, LocationInfo] = {}

def resolve_location(place: str) -> LocationInfo:
    key = place.strip()
    if not key:
        raise ValueError("place が空です")
    if key in _location_cache:
        return _location_cache[key]

    loc = _geolocator.geocode(key, language="ja")
    if loc is None:
        raise ValueError(f"場所が見つかりません: {place!r}")

    lat = float(loc.latitude)
    lon = float(loc.longitude)
    tz = _tzf.timezone_at(lat=lat, lng=lon) or "UTC"

    info = LocationInfo(name=str(loc), lat=lat, lon=lon, tz=tz)
    _location_cache[key] = info
    return info

def local_to_utc(dt_local_str: str, tz_name: str, fmt: str = "%Y/%m/%d %H:%M") -> datetime:
    dt_naive = datetime.strptime(dt_local_str, fmt)
    dt_local = dt_naive.replace(tzinfo=ZoneInfo(tz_name))
    return dt_local.astimezone(ZoneInfo("UTC"))

@dataclass
class ChartContext:
    info: LocationInfo
    dt_utc: datetime
    jd_ut: float
    elev_m: float = 0.0  # 海抜（未入力なら0でOK）

def make_context_from_place(place: str, dt_local_str: str) -> ChartContext:
    info = resolve_location(place)
    dt_utc = local_to_utc(dt_local_str, info.tz)
    jd_ut = julian_day_ut(dt_utc)
    return ChartContext(info=info, dt_utc=dt_utc, jd_ut=jd_ut, elev_m=0.0)

# ============================================================
# 計算（Swiss Ephemeris）
# ============================================================
def angular_separation_deg(a, b):
    d = abs((a - b) % 360)
    return d if d <= 180 else 360 - d

def find_aspects(longitudes, enabled_names, aspect_dict, orb_deg):
    results = []
    names = list(enabled_names)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            a, b = longitudes.get(n1), longitudes.get(n2)
            if a is None or b is None:
                continue
            sep = angular_separation_deg(a, b)

            best = None
            for asp_name, asp_deg in aspect_dict.items():
                # 特定ペアの Opposition は除外
                if asp_deg == 180 and is_excluded_opposition_pair(n1, n2):
                    continue
                delta = abs(sep - asp_deg)
                if delta <= orb_deg:
                    if (best is None) or (delta < best[-1]):
                        best = (n1, n2, asp_name, asp_deg, delta)
            if best is not None:
                results.append(best)
    return results

def deg_to_degmin(x: float):
    x = abs(x)
    d = int(x)
    m = int(round((x - d) * 60))
    if m == 60:
        d += 1
        m = 0
    return d, m

def fmt_orb(delta_deg: float) -> str:
    d, m = deg_to_degmin(delta_deg)
    return f"{d}°{m:02d}'"

def fmt_aspect_name(name: str) -> str:
    return f"{ASPECT_SYMBOLS.get(name, '?')} {name}"

def compute_planet_longitudes(jd_ut: float, body_map: dict[str, int], flags: int) -> dict[str, float | None]:
    """
    戻り値の lon は黄経（tropical, true ecliptic）を想定。
    失敗した天体は None にして、全体は落とさない。
    """
    out: dict[str, float | None] = {}
    for name, body_id in body_map.items():
        if body_id is None:
            out[name] = None
            print(f"[warn] {name}: SwissEph constant not available.")
            continue
        try:
            xx, retflag = swe.calc_ut(jd_ut, body_id, flags)
            lon = float(xx[0])
            out[name] = norm360(lon)
        except Exception as e:
            out[name] = None
            print(f"[warn] {name}: calc_ut failed ({e}).")
    return out

def compute_retrograde_flags(jd_ut: float, body_map: dict[str, int],
                             flags: int, delta_days: float = 1.0) -> dict[str, bool | None]:
    """
    前後比較による逆行判定
    - True  : 逆行
    - False : 順行
    - None  : 判定不能
    """
    out: dict[str, bool | None] = {}

    jd_next = jd_ut + delta_days

    for name, body_id in body_map.items():
        if body_id is None:
            out[name] = None
            continue

        try:
            xx0, _ = swe.calc_ut(jd_ut, body_id, flags)
            xx1, _ = swe.calc_ut(jd_next, body_id, flags)

            lon0 = norm360(float(xx0[0]))
            lon1 = norm360(float(xx1[0]))

            # 360度跨ぎ対策
            diff = (lon1 - lon0 + 540.0) % 360.0 - 180.0

            out[name] = diff < 0.0

        except Exception as e:
            print(f"[warn] {name}: retrograde check failed ({e}).")
            out[name] = None

    return out

def compute_houses_and_axes(jd_ut, lat, lon, hsys="P"):
    if isinstance(hsys, str):
        hsys = hsys.encode("ascii")

    result = swe.houses_ex(jd_ut, lat, lon, hsys)

    cusps = result[0]
    ascmc = result[1]

    asc = norm360(float(ascmc[0]))
    mc  = norm360(float(ascmc[1]))

    axes = {
        "ASC": asc,
        "DSC": norm360(asc + 180.0),
        "MC": mc,
        "IC": norm360(mc + 180.0),
        "ARMC": norm360(float(ascmc[2])),
        "Vertex": norm360(float(ascmc[3])),
    }

    # cusps は「1始まり・13要素」
    house_cusps = [norm360(float(c)) for c in cusps]

    return axes, house_cusps



def sun_altitude_deg(jd_ut: float, lon: float, lat: float, elev_m: float, flags_pos: int) -> float | None:
    """
    swe.azalt で太陽高度を出す。
    うまくいかない環境もあるので、失敗したら None を返す（=昼夜判定しない）
    """
    try:
        # equatorial coordinates for Sun
        flags_equ = (flags_pos | swe.FLG_EQUATORIAL)
        xx, _ = swe.calc_ut(jd_ut, swe.SUN, flags_equ)
        ra, dec, dist = float(xx[0]), float(xx[1]), float(xx[2])

        calc_flag = getattr(swe, "EQU2HOR", 1)  # SE_EQU2HOR
        geopos = [lon, lat, elev_m]
        xaz = swe.azalt(jd_ut, calc_flag, geopos, 0.0, 10.0, [ra, dec, dist])
        # xaz = [azimuth, true_altitude, apparent_altitude]
        return float(xaz[1])
    except Exception as e:
        print(f"[warn] Sun altitude calc failed ({e}). Day/night fallback: assume DAY.")
        return None

def part_of_fortune_lon(asc: float, sun: float, moon: float, is_day: bool) -> float:
    if is_day:
        return norm360(asc + moon - sun)
    return norm360(asc + sun - moon)

def compute_extra_points(jd_ut: float, lon: float, lat: float, elev_m: float,
                         axes: dict, longitudes: dict, flags_pos: int) -> dict:
    """
    - North Node: mean node
    - Lilith: mean lunar apogee (Black Moon)
    - Selena: 反対側（perigee扱い）として apogee + 180°
    - Vertex: houses_ex で取ったものを採用
    - PoF: 日中/夜で切替
    """
    extras: dict[str, float | bool] = {}

    # mean node
    try:
        xx, _ = swe.calc_ut(jd_ut, swe.MEAN_NODE, flags_pos)
        nn = norm360(float(xx[0]))
    except Exception as e:
        print(f"[warn] Mean node failed ({e}).")
        nn = None

    # mean apogee (Lilith)
    try:
        xx, _ = swe.calc_ut(jd_ut, swe.MEAN_APOG, flags_pos)
        apog = norm360(float(xx[0]))
    except Exception as e:
        print(f"[warn] Mean apogee failed ({e}).")
        apog = None

    # day/night
    alt = sun_altitude_deg(jd_ut, lon, lat, elev_m, flags_pos)
    is_day = True if (alt is None) else (alt > 0.0)

    # PoF needs Sun/Moon
    sun = longitudes.get("Sun")
    moon = longitudes.get("Moon")
    asc = axes["ASC"]

    pof = None
    if sun is not None and moon is not None:
        pof = part_of_fortune_lon(asc=asc, sun=sun, moon=moon, is_day=is_day)

    # Selena (perigee扱い)
    selena = None
    if apog is not None:
        selena = norm360(apog + 180.0)

    extras["NorthNode"] = nn
    extras["SouthNode"] = (norm360(nn + 180.0) if nn is not None else None)
    extras["Lilith"] = apog
    extras["Selena"] = selena
    extras["PoF"] = pof
    extras["Vertex"] = axes.get("Vertex")
    extras["IsDayChart"] = is_day
    return extras

# ============================================================
# 描画ユーティリティ
# ============================================================
def deg_to_rad(deg):
    return math.radians(deg % 360)

def draw_circle(ax, r=1.0, lw=1.0, color="black"):
    angles = [math.radians(a) for a in range(0, 361)]
    rs = [r] * len(angles)
    ax.plot(angles, rs, linewidth=lw, color=color)

def draw_radial_line(ax, deg, r0=0.0, r1=1.0, lw=1.0, color="black"):
    a = deg_to_rad(deg)
    ax.plot([a, a], [r0, r1], linewidth=lw, color=color)

# ============================================================
# GUI本体
# ============================================================
class AstroApp:
    def __init__(self, root: tk.Tk, ctx: ChartContext):
        self.root = root
        self.root.title("Astro Chart (pyswisseph)")

        self.ctx = ctx

        # ----- 計算対象 -----
        # Swiss Ephemeris の body id
        # 小惑星は ephe ファイルが無いと失敗することがあります（その場合はNoneにして続行）
        self.body_ids = {
            "Sun": swe.SUN,
            "Moon": swe.MOON,
            "Mercury": swe.MERCURY,
            "Venus": swe.VENUS,
            "Mars": swe.MARS,
            "Jupiter": swe.JUPITER,
            "Saturn": swe.SATURN,
            "Uranus": swe.URANUS,
            "Neptune": swe.NEPTUNE,
            "Pluto": swe.PLUTO,

            # --- asteroids ---
            "Ceres": getattr(swe, "CERES", 15),
            "Pallas": getattr(swe, "PALLAS", 16),
            "Juno": getattr(swe, "JUNO", 17),
            "Vesta": getattr(swe, "VESTA", 18),
            "Chiron": getattr(swe, "CHIRON", 15_000),  # fallback (実際は swe.CHIRON があるはず)
            "Pholus": getattr(swe, "PHOLUS", 15_001),
        }

        # 計算フラグ：SWIEPH（高精度）を基本に。失敗時は None で続行。
        self.flags_pos = swe.FLG_SWIEPH

        # ---------------------
        # レイアウト
        # ---------------------
        self.left = ttk.Frame(root, padding=10)
        self.left.pack(side="left", fill="y")

        self.right = ttk.Frame(root, padding=10)
        self.right.pack(side="right", fill="both", expand=True)

        # ---------------------
        # 表示トグル（惑星）
        # ---------------------
        ttk.Label(self.left, text="表示する天体", font=("", 11, "bold")).pack(anchor="w", pady=(0, 6))

        self.body_vars = {}
        for name in self.body_ids.keys():
            var = tk.BooleanVar(value=True)
            self.body_vars[name] = var
            ttk.Checkbutton(self.left, text=name, variable=var, command=self.redraw).pack(anchor="w")

        ttk.Separator(self.left).pack(fill="x", pady=8)

        # ---------------------
        # 表示トグル（計算点 / ノードなど）
        # ---------------------
        ttk.Label(self.left, text="表示する計算点", font=("", 11, "bold")).pack(anchor="w", pady=(0, 6))

        self.extra_point_names = ["NorthNode", "SouthNode", "PoF", "Lilith", "Selena", "Vertex"]
        self.extra_vars = {}
        for name in self.extra_point_names:
            var = tk.BooleanVar(value=True)
            self.extra_vars[name] = var
            ttk.Checkbutton(self.left, text=name, variable=var, command=self.redraw).pack(anchor="w")

        ttk.Separator(self.left).pack(fill="x", pady=8)

        # ---------------------
        # 表示トグル（グリッド）
        # ---------------------
        ttk.Label(self.left, text="グリッド表示", font=("", 11, "bold")).pack(anchor="w", pady=(0, 6))
        self.show_zodiac = tk.BooleanVar(value=True)
        self.show_houses = tk.BooleanVar(value=True)
        self.show_axes = tk.BooleanVar(value=True)

        ttk.Checkbutton(self.left, text="星座(12分割)", variable=self.show_zodiac, command=self.redraw).pack(anchor="w")
        ttk.Checkbutton(self.left, text="ハウス(Equal)", variable=self.show_houses, command=self.redraw).pack(anchor="w")
        ttk.Checkbutton(self.left, text="軸(AC/DC/MC/IC)", variable=self.show_axes, command=self.redraw).pack(anchor="w")

        # ---------------------
        # アスペクト表示（メジャー + マイナー）
        # ---------------------
        ttk.Separator(self.left).pack(fill="x", pady=8)
        ttk.Label(self.left, text="アスペクト", font=("", 11, "bold")).pack(anchor="w", pady=(0, 6))

        self.show_aspects = tk.BooleanVar(value=True)
        self.show_minor_aspects = tk.BooleanVar(value=False)
        ttk.Checkbutton(self.left, text="アスペクト(メジャー)", variable=self.show_aspects, command=self.redraw).pack(anchor="w")
        ttk.Checkbutton(self.left, text="マイナーも表示", variable=self.show_minor_aspects, command=self.redraw).pack(anchor="w")

        ttk.Label(self.left, text="オーブ(メジャー)", font=("", 10, "bold")).pack(anchor="w", pady=(8, 0))
        self.major_orb_var = tk.DoubleVar(value=DEFAULT_MAJOR_ORB_DEG)
        ttk.Spinbox(self.left, from_=0.5, to=12.0, increment=0.5, textvariable=self.major_orb_var,
                    width=6, command=self.redraw).pack(anchor="w")

        ttk.Label(self.left, text="オーブ(マイナー)", font=("", 10, "bold")).pack(anchor="w", pady=(6, 0))
        self.minor_orb_var = tk.DoubleVar(value=DEFAULT_MINOR_ORB_DEG)
        ttk.Spinbox(self.left, from_=0.5, to=8.0, increment=0.5, textvariable=self.minor_orb_var,
                    width=6, command=self.redraw).pack(anchor="w")

        # 日時表示
        ttk.Separator(self.left).pack(fill="x", pady=8)
        ttk.Label(self.left, text="計算日時 (UTC)", font=("", 10, "bold")).pack(anchor="w")
        self.date_label = ttk.Label(self.left, text=self.ctx.dt_utc.isoformat(timespec="seconds"))
        self.date_label.pack(anchor="w", pady=(2, 6))

        ttk.Label(self.left, text="SWE eph path", font=("", 9, "bold")).pack(anchor="w", pady=(4, 0))
        ttk.Label(self.left, text=SWE_EPH_PATH, wraplength=220).pack(anchor="w", pady=(2, 6))

        ttk.Button(self.left, text="再計算して再描画", command=self.recalc_and_redraw).pack(anchor="w")

        # ---------------------
        # matplotlib埋め込み
        # ---------------------
        self.fig, self.ax = plt.subplots(subplot_kw={'projection': 'polar'})
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        self.recalc_and_redraw()

    def recalc_and_redraw(self):
        jd_ut = self.ctx.jd_ut
        lon = self.ctx.info.lon
        lat = self.ctx.info.lat
        self.extras = {}

        jd_ut = self.ctx.jd_ut
        lon = self.ctx.info.lon
        lat = self.ctx.info.lat

        # ② 天体
        self.longitudes = compute_planet_longitudes(
            jd_ut, self.body_ids, self.flags_pos
        )

        # ③ ハウス＆軸
        hsys_code, hsys_label = choose_house_system(lat)
        self.axes, self.house_cusps = compute_houses_and_axes(
            jd_ut, lat, lon, hsys=hsys_code
        )
        self.house_system_label = hsys_label

        # ④ extras
        self.extras = compute_extra_points(
            jd_ut=jd_ut,
            lon=lon,
            lat=lat,
            elev_m=self.ctx.elev_m,
            axes=self.axes,
            longitudes=self.longitudes,
            flags_pos=self.flags_pos
        )

        # 天体経度
        self.longitudes = compute_planet_longitudes(jd_ut, self.body_ids, self.flags_pos)
        
        # 天体逆行フラグ
        self.retrogrades = compute_retrograde_flags(
        jd_ut, self.body_ids, self.flags_pos
        )

        # ハウス方式の自動選択
        hsys_code, hsys_label = choose_house_system(lat)

        # ハウス＆軸
        self.axes, self.house_cusps = compute_houses_and_axes(
            jd_ut, lat, lon, hsys=hsys_code
        )

        # extras（最終的な longitudes / axes に揃えて再計算）
        self.extras = compute_extra_points(
            jd_ut=jd_ut,
            lon=lon,
            lat=lat,
            elev_m=self.ctx.elev_m,
            axes=self.axes,
            longitudes=self.longitudes,
            flags_pos=self.flags_pos
        )

        # 表示用に保持
        self.house_system_label = hsys_label
        
        print("\n=== House System ===")
        print(self.house_system_label)


        # =====================
        # コンソール出力
        # =====================
        print("\n=== Location / Time ===")
        print(self.ctx.info)
        print("UTC:", self.ctx.dt_utc.isoformat())
        print("JD(UT):", jd_ut)

        print("\n=== Planet positions (sign / house / longitude) ===")
        ac = self.axes["ASC"]
        for name in self.body_ids.keys():
            lonv = self.longitudes.get(name)
            retro = self.retrogrades.get(name)

            rx = " ℞" if retro else ""

            if lonv is None:
                print(f"{name:8s} | (missing){rx:3s}              | House -- | Lon   --")
                continue

            sign_str = deg_to_sign_string(lonv)
            house = lon_to_house_equal(lonv, ac)

            print(f"{name:8s} | {sign_str:25s}{rx:3s} | House {house:2d} | Lon {lonv:7.2f}°")

        print("\n=== Axes (ecliptic longitude) ===")
        for k in ["ASC", "DSC", "MC", "IC"]:
            lonv = self.axes[k]
            print(f"{k:2s}       | {deg_to_sign_string(lonv):25s} | Lon {lonv:7.2f}°")
        print(f"Vertex  | {deg_to_sign_string(self.axes['Vertex']):25s} | Lon {self.axes['Vertex']:7.2f}°")

        print("\n=== Extras ===")
        for k in ["NorthNode", "SouthNode", "PoF", "Lilith", "Selena", "Vertex"]:
            lonv = self.extras.get(k)
            if lonv is None:
                print(f"{k:9s} | (missing)                 | House -- | Lon   --")
                continue
            print(f"{k:9s} | {deg_to_sign_string(lonv):25s} | House {lon_to_house_equal(lonv, ac):2d} | Lon {lonv:7.2f}°")
        print("Day chart?", self.extras.get("IsDayChart"))

        self.date_label.config(text=self.ctx.dt_utc.isoformat(timespec="seconds"))

        # =====================
        # コンソール出力：アスペクト
        # =====================
        # =====================
        # アスペクト/パターン判定用：惑星 + 軸 + 計算点
        # =====================
        self.all_longitudes = dict(self.longitudes)

        # Axes / Angles
        for k in ["ASC", "DSC", "MC", "IC", "Vertex"]:
            self.all_longitudes[k] = self.axes.get(k)

        # Extras（IsDayChart は角度ではないので除外）
        for k in self.extra_point_names:
            v = self.extras.get(k)
            self.all_longitudes[k] = (float(v) % 360.0) if isinstance(v, (int, float)) and v is not None else None

        enabled = []

        # planets（チェックボックスに従う）
        for n in self.body_ids.keys():
            if self.body_vars[n].get() and self.all_longitudes.get(n) is not None:
                enabled.append(n)

        # axes（常に含める）
        for n in ["ASC", "DSC", "MC", "IC", "Vertex"]:
            if self.all_longitudes.get(n) is not None:
                enabled.append(n)

        # extras（チェックボックスに従う）
        for n in self.extra_point_names:
            if self.extra_vars.get(n).get() and self.all_longitudes.get(n) is not None:
                enabled.append(n)

        # 重複除去（順序維持）
        enabled = list(dict.fromkeys(enabled))

        major_orb = float(self.major_orb_var.get())
        minor_orb = float(self.minor_orb_var.get())

        major_hits = find_aspects(self.all_longitudes, enabled, MAJOR_ASPECTS, orb_deg=major_orb)
        minor_hits = find_aspects(self.all_longitudes, enabled, MINOR_ASPECTS, orb_deg=minor_orb)
        major_hits.sort(key=lambda t: t[-1])
        minor_hits.sort(key=lambda t: t[-1])

        print("\n=== Aspects (MAJOR) ===")
        print(f"(orb <= {major_orb}°) targets={len(enabled)} bodies hits={len(major_hits)}")
        if not major_hits:
            print("  (none)")
        else:
            for n1, n2, asp_name, asp_deg, delta in major_hits:
                print(f"  {n1:8s} {fmt_aspect_name(asp_name):18s} {n2:8s} | exact {asp_deg:3d}° | orb {fmt_orb(delta)}")

        print("\n=== Aspects (MINOR) ===")
        print(f"(orb <= {minor_orb}°) targets={len(enabled)} bodies hits={len(minor_hits)}")
        if not minor_hits:
            print("  (none)")
        else:
            for n1, n2, asp_name, asp_deg, delta in minor_hits:
                print(f"  {n1:8s} {fmt_aspect_name(asp_name):18s} {n2:8s} | exact {asp_deg:3d}° | orb {fmt_orb(delta)}")
        
        # =====================
        # コンソール出力：チャートパターン（追加）
        # =====================
        patterns = detect_chart_patterns(self.all_longitudes, enabled, major_orb, minor_orb)
        print_chart_patterns(patterns)

        self.redraw()

    def redraw(self):
        self.ax.clear()

        # polar設定
        self.ax.set_theta_zero_location("E")
        self.ax.set_theta_direction(-1)
        self.ax.set_rlim(0, 1.20)
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])

        # 外円
        draw_circle(self.ax, r=1.00, lw=1.0, color="black")

        # 星座12分割（絶対座標）
        if self.show_zodiac.get():
            for i in range(12):
                deg = i * 30
                draw_radial_line(self.ax, deg, r0=0.0, r1=1.0, lw=0.6, color="gray")
                label_deg = deg + 15
                self.ax.text(deg_to_rad(label_deg), 1.10, SIGN_SYMBOLS[i],
                             ha="center", va="center", fontsize=12)

        # ハウス（Equal, AC基準）
        if self.show_houses.get():
            for i, cusp_deg in enumerate(self.house_cusps):
                draw_radial_line(self.ax, cusp_deg, r0=0.0, r1=1.0, lw=1.0, color="dimgray")
                mid_deg = (cusp_deg + 15) % 360
                self.ax.text(deg_to_rad(mid_deg), 0.78, str(i + 1),
                             ha="center", va="center", fontsize=9)

        # 軸（AC/DC/MC/IC）
        if self.show_axes.get():
            for key in ["ASC", "DSC", "MC", "IC"]:
                draw_radial_line(self.ax, self.axes[key], r0=0.0, r1=1.0, lw=2.2, color="black")
            for key in ["ASC", "DSC", "MC", "IC"]:
                self.ax.text(deg_to_rad(self.axes[key]), 1.16, key,
                             ha="center", va="center", fontsize=9)

        # アスペクト線
        if self.show_aspects.get():
            enabled = []

            # planets（チェックボックス）
            for n in self.body_ids.keys():
                if self.body_vars[n].get() and self.all_longitudes.get(n) is not None:
                    enabled.append(n)

            # axes（常に）
            for n in ["ASC", "DSC", "MC", "IC", "Vertex"]:
                if self.all_longitudes.get(n) is not None:
                    enabled.append(n)

            # extras（チェックボックス）
            for n in self.extra_point_names:
                if self.extra_vars.get(n).get() and self.all_longitudes.get(n) is not None:
                    enabled.append(n)

            enabled = list(dict.fromkeys(enabled))
            r_aspect = 0.65

            major_orb = float(self.major_orb_var.get())
            major_hits = find_aspects(self.all_longitudes, enabled, MAJOR_ASPECTS, orb_deg=major_orb)
            for n1, n2, asp_name, asp_deg, delta in major_hits:
                a1 = deg_to_rad(self.all_longitudes[n1])
                a2 = deg_to_rad(self.all_longitudes[n2])
                color = ASPECT_COLORS.get(asp_name, "gray")
                self.ax.plot([a1, a2], [r_aspect, r_aspect],
                             linewidth=1.3, color=color, alpha=0.9, zorder=2)

            if self.show_minor_aspects.get():
                minor_orb = float(self.minor_orb_var.get())
                minor_hits = find_aspects(self.all_longitudes, enabled, MINOR_ASPECTS, orb_deg=minor_orb)
                for n1, n2, asp_name, asp_deg, delta in minor_hits:
                    a1 = deg_to_rad(self.all_longitudes[n1])
                    a2 = deg_to_rad(self.all_longitudes[n2])
                    color = ASPECT_COLORS.get(asp_name, "gray")
                    self.ax.plot([a1, a2], [r_aspect, r_aspect],
                                 linewidth=0.9, color=color, alpha=0.55, zorder=1)

        # 天体プロット
        for name in self.body_ids.keys():
            if not self.body_vars[name].get():
                continue
            lonv = self.longitudes.get(name)
            if lonv is None:
                continue
            ang = deg_to_rad(lonv)
            self.ax.scatter([ang], [1.00], s=40, zorder=5)
            self.ax.text(ang, 1.05, name, ha="center", va="center", fontsize=9)

        # 計算点プロット（ノード/リリス/セレナ/PoF/Vertex）
        for name in self.extra_point_names:
            if not self.extra_vars.get(name).get():
                continue
            lonv = self.all_longitudes.get(name)
            if lonv is None:
                continue
            ang = deg_to_rad(lonv)
            self.ax.scatter([ang], [0.92], s=28, zorder=5)
            self.ax.text(ang, 0.97, name, ha="center", va="center", fontsize=8)

        self.canvas.draw()

# ============================================================
# エントリポイント
# ============================================================
if __name__ == "__main__":
    print("=== astromaker (pyswisseph) ===")
    print("Tip: 小惑星を正確に出すには ephemeris ファイル（seas_18.se1 など）を SWE_EPH_PATH に置くと安定します。")
    print("SWE_EPH_PATH =", SWE_EPH_PATH)

    place = input("場所を入力 (例: Tokyo, Japan / 名古屋): ").strip() or "Tokyo, Japan"
    dt_local = input("ローカル時刻を入力 (例: 2001/01/01 00:00): ").strip() or "2009/10/30 11:32"

    try:
        ctx = make_context_from_place(place, dt_local)
        print("\n=== Location Applied ===")
        print(ctx.info)
        print("UTC:", ctx.dt_utc.isoformat())
        print("JD(UT):", ctx.jd_ut)
    except Exception as e:
        print("\n場所/時刻の解決に失敗:", e)
        print("デフォルト（Tokyo, Japan / 2001/01/01 00:00）で続行します。")
        ctx = make_context_from_place("Tokyo, Japan", "2001/01/01 00:00")

    root = tk.Tk()
    app = AstroApp(root, ctx)
    root.mainloop()
