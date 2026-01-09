import tkinter as tk
from tkinter import ttk

import ephem
import math

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# =====================
# 定数（星座）
# =====================
SIGNS = ["♈︎","♉︎","♊︎","♋︎","♌︎","♍︎","♎︎","♏︎","♐︎","♑︎","♒︎","♓︎"]

# =====================
# アスペクト設定（基礎5つ）
# =====================
ASPECTS = {
    "Conjunction": 0,
    "Sextile": 60,
    "Square": 90,
    "Trine": 120,
    "Opposition": 180,
}

# 種類ごとの色（見やすさ優先で固定）
ASPECT_COLORS = {
    "Conjunction": "gray",
    "Sextile": "dodgerblue",
    "Square": "red",
    "Trine": "limegreen",
    "Opposition": "orange",
}

# 角度の許容誤差（オーブ）
DEFAULT_ORB_DEG = 6.0

# 計算エリア
def compute_longitudes(obs, planets_dict):
    """惑星の黄経（度, 0-360）を計算"""
    longitudes = {}
    for name, body in planets_dict.items():
        body.compute(obs)
        lon = math.degrees(body.hlong) % 360
        longitudes[name] = lon
    return longitudes

def angular_separation_deg(a, b):
    """角度差（0〜180度）"""
    d = abs((a - b) % 360)
    return d if d <= 180 else 360 - d


def find_aspects(longitudes, enabled_names, orb_deg=DEFAULT_ORB_DEG):
    """
    longitudes: dict[name] = deg
    enabled_names: 表示対象の天体名リスト
    戻り値: list of (name1, name2, aspect_name, exact_deg, delta)
    """
    results = []
    names = list(enabled_names)

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            a, b = longitudes[n1], longitudes[n2]
            sep = angular_separation_deg(a, b)

            # 基礎5つのどれに近いか
            for asp_name, asp_deg in ASPECTS.items():
                delta = abs(sep - asp_deg)
                if delta <= orb_deg:
                    results.append((n1, n2, asp_name, asp_deg, delta))
                    break

    return results

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

        # オーブ
        ttk.Label(self.left, text="オーブ(度)", font=("", 10, "bold")).pack(anchor="w", pady=(8, 0))
        self.orb_var = tk.DoubleVar(value=DEFAULT_ORB_DEG)
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

        # デバッグ出力（必要なら消してOK）
        print("=== Longitudes (deg) ===")
        for k, v in self.longitudes.items():
            print(f"{k:8s}: {v:7.2f}°")
        print("=== Axes (deg) ===", {k: round(v, 2) for k, v in self.axes.items()})
        print("=== House cusps (deg) ===", [round(x, 2) for x in self.house_cusps])

        self.date_label.config(text=str(self.obs.date))
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
                    SIGNS[i],
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
        if self.show_aspects.get():
            # いま表示ONの天体だけを対象にする（チェックと一致）
            enabled = [name for name in self.longitudes.keys() if self.body_vars[name].get()]
            orb = float(self.orb_var.get()) if hasattr(self, "orb_var") else DEFAULT_ORB_DEG

            aspects = find_aspects(self.longitudes, enabled, orb_deg=orb)

            # 線は内側の半径で結ぶ（外周だとゴチャる）
            r_aspect = 0.80

            for n1, n2, asp_name, asp_deg, delta in aspects:
                a1 = deg_to_rad(self.longitudes[n1])
                a2 = deg_to_rad(self.longitudes[n2])

                color = ASPECT_COLORS.get(asp_name, "gray")
                lw = 1.2

                # 2点を直線で結ぶ（polar上でもOK）
                self.ax.plot([a1, a2], [r_aspect, r_aspect], linewidth=lw, color=color, alpha=0.85, zorder=2)


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
