import tkinter as tk
from tkinter import ttk

import math
import os

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

        # 天体経度
        self.longitudes = compute_planet_longitudes(jd_ut, self.body_ids, self.flags_pos)

        # ハウス＆軸（Equal）
        self.axes, self.house_cusps = compute_houses_and_axes(jd_ut, lat, lon, hsys="P")

        # extras（ノード/PoF/Vertex/Lilith/Selena）
        self.extras = compute_extra_points(
            jd_ut=jd_ut, lon=lon, lat=lat, elev_m=self.ctx.elev_m,
            axes=self.axes, longitudes=self.longitudes, flags_pos=self.flags_pos
        )

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
            if lonv is None:
                print(f"{name:8s} | (missing)                 | House -- | Lon   --")
                continue
            sign_str = deg_to_sign_string(lonv)
            house = lon_to_house_equal(lonv, ac)
            print(f"{name:8s} | {sign_str:25s} | House {house:2d} | Lon {lonv:7.2f}°")

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
        enabled = [n for n in self.body_ids.keys() if self.body_vars[n].get() and self.longitudes.get(n) is not None]
        major_orb = float(self.major_orb_var.get())
        minor_orb = float(self.minor_orb_var.get())

        major_hits = find_aspects(self.longitudes, enabled, MAJOR_ASPECTS, orb_deg=major_orb)
        minor_hits = find_aspects(self.longitudes, enabled, MINOR_ASPECTS, orb_deg=minor_orb)
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
            enabled = [n for n in self.body_ids.keys() if self.body_vars[n].get() and self.longitudes.get(n) is not None]
            r_aspect = 0.65

            major_orb = float(self.major_orb_var.get())
            major_hits = find_aspects(self.longitudes, enabled, MAJOR_ASPECTS, orb_deg=major_orb)
            for n1, n2, asp_name, asp_deg, delta in major_hits:
                a1 = deg_to_rad(self.longitudes[n1])
                a2 = deg_to_rad(self.longitudes[n2])
                color = ASPECT_COLORS.get(asp_name, "gray")
                self.ax.plot([a1, a2], [r_aspect, r_aspect],
                             linewidth=1.3, color=color, alpha=0.9, zorder=2)

            if self.show_minor_aspects.get():
                minor_orb = float(self.minor_orb_var.get())
                minor_hits = find_aspects(self.longitudes, enabled, MINOR_ASPECTS, orb_deg=minor_orb)
                for n1, n2, asp_name, asp_deg, delta in minor_hits:
                    a1 = deg_to_rad(self.longitudes[n1])
                    a2 = deg_to_rad(self.longitudes[n2])
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
