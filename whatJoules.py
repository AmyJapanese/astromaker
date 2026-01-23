# タイピング仕事量計算ツール
# 前提：理想キーボード（1打鍵 = 0.005 J）

ENERGY_PER_KEY = 0.005  # J
PERSON_WEIGHT = 50      # kg
GRAVITY = 9.8           # m/s^2
CAT_WEIGHT = 10 #猫の重さ

# 入力
keystrokes = int(input("今日のタイピング回数を入力してください: "))

# 計算
energy_joule = keystrokes * ENERGY_PER_KEY
height_m = energy_joule / (PERSON_WEIGHT * GRAVITY)
height_m_cat = energy_joule / (CAT_WEIGHT * GRAVITY)

# 出力
print("\n--- 今日のタイピング成果 ---")
print(f"打鍵数: {keystrokes:,} 回")
print(f"発生した仕事量: {energy_joule:.1f} J")
print(f"{PERSON_WEIGHT}kgの人を持ち上げた高さ: {height_m:.3f} m")
print(f"（約 {height_m*100:.1f} cm）")

print(f"{CAT_WEIGHT}kgの猫を持ち上げた高さ: {height_m_cat:.3f} m")
print(f"（約 {height_m_cat*100:.1f} cm）")