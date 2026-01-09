import tkinter as tk
from tkinter import ttk

import math
import ephem

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# =====================
# 定数（星座）
# =====================
SIGN_NAMES = ["Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
              "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius", "Pisces"]

SIGN_SYMBOLS = ["♈︎","♉︎","♊︎","♋︎","♌︎","♍︎","♎︎","♏︎","♐︎","♑︎","♒︎","♓︎"]

def deg_to_sign_index(lon_deg: float) -> int:
    """黄経(0-360) → 星座index(0-11)"""
    return int((lon_deg % 360) // 30)

def deg_to_sign_pos(lon_deg: float):
    """
    黄経(0-360) → (星座index, 星座内度数(0-30), 分)
    分まで欲しいならこの形が便利
    """
    lon = lon_deg % 360
    sign_idx = deg_to_sign_index(lon)
    in_sign = lon - sign_idx * 30
    deg_int = int(in_sign)
    minute = int(round((in_sign - deg_int) * 60))
    # 60分になったときの繰り上げケア（丸め誤差対策）
    if minute == 60:
        deg_int += 1
        minute = 0
        if deg_int == 30:
            deg_int = 0
            sign_idx = (sign_idx + 1) % 12
    return sign_idx, deg_int, minute

def deg_to_sign_string(lon_deg: float) -> str:
    """例: ♑︎ 12°34' (Capricorn) みたいな表示文字列"""
    s, d, m = deg_to_sign_pos(lon_deg)
    return f"{SIGN_SYMBOLS[s]} {d:02d}°{m:02d}' ({SIGN_NAMES[s]})"

def lon_to_house_equal(lon_deg: float, ac_deg: float) -> int:
    """
    Equal House専用：ACを1ハウス起点として30°刻み
    """
    rel = (lon_deg - ac_deg) % 360
    return int(rel // 30) + 1


# =====================
# アスペクト（メジャー/マイナー）
# =====================
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

# 色（見やすさ優先：好きに変えてOK）
ASPECT_COLORS = {
    # Major
    "Conjunction": "gray",
    "Sextile": "dodgerblue",
    "Square": "red",
    "Trine": "limegreen",
    "Opposition": "orange",

    # Minor（目立ちすぎない色にしておくと崩れにくい）
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

# オーブ（メジャー/マイナーで分けるのが実用的）
DEFAULT_MAJOR_ORB_DEG = 6.0
DEFAULT_MINOR_ORB_DEG = 2.0

# 計算エリア
def compute_longitudes(obs, planets_dict):
    """惑星の黄経（地球中心・黄道座標 / 度, 0-360）を計算"""
    longitudes = {}
    for name, body in planets_dict.items():
        body.compute(obs)
        ecl = ephem.Ecliptic(body)                 # ← これが地球中心の黄道座標に変換してくれる
        lon = math.degrees(float(ecl.lon)) % 360   # lonはラジアンなので必ずdegrees
        longitudes[name] = lon
    return longitudes

def angular_separation_deg(a, b):
    """角度差（0〜180度）"""
    d = abs((a - b) % 360)
    return d if d <= 180 else 360 - d


def find_aspects(longitudes, enabled_names, aspect_dict, orb_deg):
    """
    longitudes: dict[name] = deg
    enabled_names: 表示対象の天体名リスト
    aspect_dict: {"Square": 90, ...} みたいな辞書
    orb_deg: 許容誤差（度）
    戻り値: list of (name1, name2, aspect_name, exact_deg, delta)
    """
    results = []
    names = list(enabled_names)

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            a, b = longitudes[n1], longitudes[n2]
            sep = angular_separation_deg(a, b)

            # 近い順に当てる（同じorb内で複数に引っかかるのを避ける）
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
    """度(小数) → (度, 分)"""
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

def calc_axes_degrees(obs):
    """
    AC/DC/MC/IC を度(0-360)で返す（黄経）
    obs.sidereal_time() はラジアン（＝赤経系の角度）なので、
    占星術の軸（黄経）へ変換する。
    """
    theta = float(obs.sidereal_time())     # 地方恒星時 θ [rad]
    phi = float(obs.lat)                   # 緯度 φ [rad]
    eps = math.radians(23.439291)          # 黄道傾斜角 ε [rad]

    # --- MC（黄経） ---
    # λ_MC = atan2( sinθ, cosθ * cosε )
    mc = math.degrees(math.atan2(math.sin(theta), math.cos(theta) * math.cos(eps))) % 360

    # --- AC（黄経） ---
    # よく使われる形： λ_ASC = atan2( -cosθ, sinθ*cosε + tanφ*sinε ) + 180°
    asc = math.degrees(math.atan2(
        -math.cos(theta),
        (math.sin(theta) * math.cos(eps) + math.tan(phi) * math.sin(eps))
    ))
    ac = (asc + 180.0) % 360

    return {
        "AC": ac,
        "DC": (ac + 180) % 360,
        "MC": mc,
        "IC": (mc + 180) % 360
    }



def calc_equal_house_cusps(ac_deg):
    """Equal House：1Hカスプ=AC、30°ずつの12カスプ（度）"""
    return [(ac_deg + i * 30) % 360 for i in range(12)]


# =====================
# 描画ユーティリティ
# =====================
def deg_to_rad(deg):
    return math.radians(deg % 360)


def draw_circle(ax, r=1.0, lw=1.0, color="black"):
    angles = [math.radians(a) for a in range(0, 361)]
    rs = [r] * len(angles)
    ax.plot(angles, rs, linewidth=lw, color=color)


def draw_radial_line(ax, deg, r0=0.0, r1=1.0, lw=1.0, color="black"):
    a = deg_to_rad(deg)
    ax.plot([a, a], [r0, r1], linewidth=lw, color=color)


# =====================
# GUI本体
# =====================
class AstroApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Astro Chart (Ephem + Tkinter + Matplotlib)")

        # ---------------------
        # 観測地点＆日時（UTC）
        # ---------------------
        self.obs = ephem.Observer()
        self.obs.lat = '35.6895'
        self.obs.lon = '139.6917'
        self.obs.date = '2009/10/30 02:32'  # UTC扱い（あなたの値）

        # ---------------------
        # 天体
        # ---------------------
        self.planets = {
            'Sun': ephem.Sun(),
            'Moon': ephem.Moon(),
            'Mercury': ephem.Mercury(),
            'Venus': ephem.Venus(),
            'Mars': ephem.Mars(),
            'Jupiter': ephem.Jupiter(),
            'Saturn': ephem.Saturn(),
            'Uranus': ephem.Uranus(),
            'Neptune': ephem.Neptune(),
            'Pluto': ephem.Pluto()
        }

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
        for name in self.planets.keys():
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
        
        # アスペクト表示
        self.show_aspects = tk.BooleanVar(value=True)
        ttk.Checkbutton(self.left, text="アスペクト(基礎5つ)", variable=self.show_aspects, command=self.redraw).pack(anchor="w")

        # アスペクト表示
        self.show_aspects = tk.BooleanVar(value=True)         # 線そのもののON/OFF
        self.show_minor_aspects = tk.BooleanVar(value=False)  # マイナーはデフォルトOFF

        ttk.Checkbutton(self.left, text="アスペクト(メジャー)", variable=self.show_aspects, command=self.redraw).pack(anchor="w")
        ttk.Checkbutton(self.left, text="マイナーも表示", variable=self.show_minor_aspects, command=self.redraw).pack(anchor="w")

        # オーブ（メジャー/マイナーで別にしておくと神）
        ttk.Label(self.left, text="オーブ(メジャー)", font=("", 10, "bold")).pack(anchor="w", pady=(8, 0))
        self.major_orb_var = tk.DoubleVar(value=DEFAULT_MAJOR_ORB_DEG)
        ttk.Spinbox(self.left, from_=0.5, to=12.0, increment=0.5, textvariable=self.major_orb_var, width=6, command=self.redraw).pack(anchor="w")

        ttk.Label(self.left, text="オーブ(マイナー)", font=("", 10, "bold")).pack(anchor="w", pady=(6, 0))
        self.minor_orb_var = tk.DoubleVar(value=DEFAULT_MINOR_ORB_DEG)
        ttk.Spinbox(self.left, from_=0.5, to=8.0, increment=0.5, textvariable=self.minor_orb_var, width=6, command=self.redraw).pack(anchor="w")

        # オーブ
        ttk.Label(self.left, text="オーブ(度)", font=("", 10, "bold")).pack(anchor="w", pady=(8, 0))
        self.orb_var = tk.DoubleVar(value=DEFAULT_MAJOR_ORB_DEG)
        ttk.Spinbox(self.left, from_=0.5, to=12.0, increment=0.5, textvariable=self.orb_var, width=6, command=self.redraw).pack(anchor="w")

        # 日時表示
        ttk.Label(self.left, text="計算日時 (UTC)", font=("", 10, "bold")).pack(anchor="w")
        self.date_label = ttk.Label(self.left, text=str(self.obs.date))
        self.date_label.pack(anchor="w", pady=(2, 6))

        # 再計算
        ttk.Button(self.left, text="再計算して再描画", command=self.recalc_and_redraw).pack(anchor="w")

        # ---------------------
        # matplotlib埋め込み
        # ---------------------
        self.fig, self.ax = plt.subplots(subplot_kw={'projection': 'polar'})
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # 初回計算＆描画
        self.recalc_and_redraw()

    def recalc_and_redraw(self):
        # 計算
        self.longitudes = compute_longitudes(self.obs, self.planets)
        self.axes = calc_axes_degrees(self.obs)
        self.house_cusps = calc_equal_house_cusps(self.axes["AC"])
        # =====================
        # コンソール出力：星座・度数・ハウス
        # =====================
        print("\n=== Planet positions (sign / house / longitude) ===")
        ac = self.axes["AC"]

    # 表示順を固定したいならリストで回すのが安心
        order = list(self.planets.keys())

        for name in order:
            lon = self.longitudes[name]
            sign_str = deg_to_sign_string(lon)
            house = lon_to_house_equal(lon, ac)
            print(f"{name:8s} | {sign_str:25s} | House {house:2d} | Lon {lon:7.2f}°")

        self.date_label.config(text=str(self.obs.date))
        
        # =====================
        # コンソール出力：アスペクト一覧
        # =====================
        enabled = [n for n in self.longitudes.keys() if self.body_vars[n].get()]

        major_orb = float(self.major_orb_var.get()) if hasattr(self, "major_orb_var") else DEFAULT_MAJOR_ORB_DEG
        minor_orb = float(self.minor_orb_var.get()) if hasattr(self, "minor_orb_var") else DEFAULT_MINOR_ORB_DEG

        major_hits = find_aspects(self.longitudes, enabled, MAJOR_ASPECTS, orb_deg=major_orb)
        # オーブが小さい（＝タイト）順に並べる
        major_hits.sort(key=lambda t: t[-1])

        print("\n=== Aspects (MAJOR) ===")
        print(f"(orb <= {major_orb}°)  targets={len(enabled)} bodies")
        if not major_hits:
            print("  (none)")
        else:
            for n1, n2, asp_name, asp_deg, delta in major_hits:
                symname = fmt_aspect_name(asp_name)
                print(f"  {n1:8s} {symname:18s} {n2:8s} | exact {asp_deg:3d}° | orb {fmt_orb(delta)}")

        # マイナーは「チェックONの時だけ出す」設計にするならこれ
        if getattr(self, "show_minor_aspects", tk.BooleanVar(value=False)).get():
            minor_hits = find_aspects(self.longitudes, enabled, MINOR_ASPECTS, orb_deg=minor_orb)
            minor_hits.sort(key=lambda t: t[-1])

            print("\n=== Aspects (MINOR) ===")
            print(f"(orb <= {minor_orb}°)  targets={len(enabled)} bodies")
            if not minor_hits:
                print("  (none)")
            else:
                for n1, n2, asp_name, asp_deg, delta in minor_hits:
                    symname = fmt_aspect_name(asp_name)
                    print(f"  {n1:8s} {symname:18s} {n2:8s} | exact {asp_deg:3d}° | orb {fmt_orb(delta)}")

        self.redraw()

    def redraw(self):
        # ここが「毎回書き直すキャンバス」：必ずclearする
        self.ax.clear()

        # polar設定
        self.ax.set_theta_zero_location("E")   # 0°を右
        self.ax.set_theta_direction(-1)        # 時計回り
        self.ax.set_rlim(0, 1.20)
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])

        # 外円
        draw_circle(self.ax, r=1.00, lw=1.0, color="black")

        # ---------------------
        # 星座12分割（絶対座標：0°,30°,60°...）
        # ---------------------
        if self.show_zodiac.get():
            for i in range(12):
                deg = i * 30
                draw_radial_line(self.ax, deg, r0=0.0, r1=1.0, lw=0.6, color="gray")

                # 星座記号（各30°の中央）
                label_deg = deg + 15
                self.ax.text(
                    deg_to_rad(label_deg),
                    1.10,
                    SIGN_SYMBOLS[i],
                    ha="center", va="center",
                    fontsize=12
                )

        # ---------------------
        # ハウス12分割（相対座標：AC基準）
        # ---------------------
        if self.show_houses.get():
            for i, cusp_deg in enumerate(self.house_cusps):
                draw_radial_line(self.ax, cusp_deg, r0=0.0, r1=1.0, lw=1.0, color="dimgray")

                # ハウス番号（各ハウス中央）
                mid_deg = (cusp_deg + 15) % 360
                self.ax.text(
                    deg_to_rad(mid_deg),
                    0.78,
                    str(i + 1),
                    ha="center", va="center",
                    fontsize=9
                )

        # ---------------------
        # 軸（AC/DC/MC/IC）太線
        # ---------------------
        if self.show_axes.get():
            # 太くする
            for key in ["AC", "DC", "MC", "IC"]:
                draw_radial_line(self.ax, self.axes[key], r0=0.0, r1=1.0, lw=2.2, color="black")

            # ラベル（外側に）
            for key in ["AC", "DC", "MC", "IC"]:
                self.ax.text(
                    deg_to_rad(self.axes[key]),
                    1.16,
                    key,
                    ha="center", va="center",
                    fontsize=9
                )
        
        # ---------------------
        # アスペクト（基礎5つ）線描画
        # ---------------------
        # ---------------------
        # アスペクト線描画（メジャー + 任意でマイナー）
        # ---------------------
        if self.show_aspects.get():
            enabled = [name for name in self.longitudes.keys() if self.body_vars[name].get()]

            r_aspect = 0.65

            # まずメジャー
            major_orb = float(self.major_orb_var.get())
            major_hits = find_aspects(self.longitudes, enabled, MAJOR_ASPECTS, orb_deg=major_orb)

            for n1, n2, asp_name, asp_deg, delta in major_hits:
                a1 = deg_to_rad(self.longitudes[n1])
                a2 = deg_to_rad(self.longitudes[n2])
                color = ASPECT_COLORS.get(asp_name, "gray")
                self.ax.plot([a1, a2], [r_aspect, r_aspect],
                            linewidth=1.3, color=color, alpha=0.85, zorder=2)

            # マイナー（チェックONの時だけ）
            if self.show_minor_aspects.get():
                minor_orb = float(self.minor_orb_var.get())
                minor_hits = find_aspects(self.longitudes, enabled, MINOR_ASPECTS, orb_deg=minor_orb)

                for n1, n2, asp_name, asp_deg, delta in minor_hits:
                    a1 = deg_to_rad(self.longitudes[n1])
                    a2 = deg_to_rad(self.longitudes[n2])
                    color = ASPECT_COLORS.get(asp_name, "gray")
                    # マイナーは少し薄く・細く
                    self.ax.plot([a1, a2], [r_aspect, r_aspect],
                                linewidth=0.9, color=color, alpha=0.55, zorder=1)

        # ---------------------
        # 惑星（絶対座標：黄経） ※scatterで線が絶対出ない
        # ---------------------
        for name, lon in self.longitudes.items():
            if not self.body_vars[name].get():
                continue

            ang = deg_to_rad(lon)

            # 点
            self.ax.scatter([ang], [1.00], s=40, zorder=5)

            # ラベル（少し外）
            self.ax.text(
                ang,
                1.05,
                name,
                ha="center", va="center",
                fontsize=9
            )

        self.canvas.draw()


if __name__ == "__main__":
    root = tk.Tk()
    app = AstroApp(root)
    root.mainloop()
