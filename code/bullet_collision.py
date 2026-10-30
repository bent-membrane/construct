# -*- coding: utf-8 -*-
"""
图 05：子弹星团碰撞（动态）
二维膜上两个"中心点 + 雾"相向而行：
  - 中心点（气体）：碰撞时受摩擦减速
  - 雾（无碰撞的投影）：直接穿透，保持匀速
t = 0 → 150 Myr
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- 单位与参数
KPC_KM = 3.0857e16      # 1 kpc = 3.0857e16 km
YR_S = 3.1558e7         # 1 yr = 3.1558e7 s

V0 = 2350.0             # km/s，每个团的初速度（相对 4700）
X0 = 117.5              # kpc，初始离中心距离
TAU = 30e6              # yr，气体减速的时间常数
T_END = 150e6           # yr

# 碰撞时刻：两个团各走 X0
t_c = X0 * KPC_KM / V0 / YR_S
print(f"碰撞时刻 t_c = {t_c/1e6:.1f} Myr")

KM_PER_MYR = V0 * YR_S * 1e6 / KPC_KM   # 1 Myr 走多少 kpc（用于输出）


def traj(t):
    """团 A（从 -X0 出发向右）；返回 (气体位置, 雾位置)，单位 kpc。"""
    xm = -X0 + V0 * t * YR_S / KPC_KM
    if t <= t_c:
        return xm, xm
    xg_c = -X0 + V0 * t_c * YR_S / KPC_KM
    dx = V0 * TAU * YR_S / KPC_KM * (1 - np.exp(-(t - t_c) / TAU))
    return xg_c + dx, xm


# ---------------------------------------------------------------- 画布
fig, ax = plt.subplots(figsize=(9, 4.8), facecolor="#05060a")
ax.set_facecolor("#05060a")
ax.set_xlim(-480, 480)
ax.set_ylim(-170, 170)
ax.set_aspect("equal")
ax.set_xlabel("膜上位置（kpc）", color="#bbbbbb")
ax.set_title("子弹星团：气体撞停，雾穿透", color="white", fontsize=13, pad=10)
ax.tick_params(colors="#8899aa")
for s in ax.spines.values():
    s.set_color("#223344")

# 粒子云（代表"雾"的投影）
rng = np.random.default_rng(2)
N_P = 1400
dA_x, dA_y = rng.normal(0, 30, N_P), rng.normal(0, 30, N_P)
dB_x, dB_y = rng.normal(0, 30, N_P), rng.normal(0, 30, N_P)

mist_A = ax.scatter([], [], s=2.2, c="#3fa9f5", alpha=0.10, linewidths=0)
mist_B = ax.scatter([], [], s=2.2, c="#3fa9f5", alpha=0.10, linewidths=0)

# 气体：外晕 + 亮核
glow_A = ax.scatter([], [], s=260, c="#ff7b00", alpha=0.18, linewidths=0)
glow_B = ax.scatter([], [], s=260, c="#ff7b00", alpha=0.18, linewidths=0)
core_A = ax.scatter([], [], s=60, c="#ffd166", alpha=0.95, linewidths=0)
core_B = ax.scatter([], [], s=60, c="#ffd166", alpha=0.95, linewidths=0)

txt_t = ax.text(0.02, 0.95, "", transform=ax.transAxes, color="white",
                fontsize=16, va="top", family="monospace")
txt_d = ax.text(0.02, 0.06, "", transform=ax.transAxes, color="#9fb3c8",
                fontsize=10, va="bottom")

# 图例
ax.scatter([], [], s=40, c="#ffd166", label="气体（中心点）")
ax.scatter([], [], s=40, c="#3fa9f5", alpha=0.5, label="雾（投影）")
ax.legend(loc="upper right", frameon=False, labelcolor="white", fontsize=9)

FRAMES = np.arange(0, T_END + 1, 5e6)


def update(t):
    xgA, xmA = traj(t)
    xgB, xmB = -xgA, -xmA

    mist_A.set_offsets(np.column_stack([xmA + dA_x, dA_y]))
    mist_B.set_offsets(np.column_stack([xmB + dB_x, -dA_y]))
    for obj, x in [(glow_A, xgA), (glow_B, xgB), (core_A, xgA), (core_B, xgB)]:
        obj.set_offsets([[x, 0.0]])

    txt_t.set_text(f"t = {t/1e6:5.0f} Myr")
    txt_d.set_text(
        f"气体间距 {abs(xgA-xgB):6.0f} kpc   |   雾间距 {abs(xmA-xmB):6.0f} kpc"
    )
    return mist_A, mist_B, glow_A, glow_B, core_A, core_B, txt_t, txt_d


ani = FuncAnimation(fig, update, frames=FRAMES, interval=120, blit=False)
p = os.path.join(OUT, "05_bullet_collision.gif")
ani.save(p, writer=PillowWriter(fps=8), dpi=110)
plt.close()
print("saved", p, f"{os.path.getsize(p)/1e6:.2f} MB")

# 关键帧数值
print("\n关键帧：")
for t in [0, 50e6, t_c, 100e6, 150e6]:
    xg, xm = traj(t)
    print(f"  t={t/1e6:5.1f} Myr  气体 x=±{abs(xg):6.1f}  雾 x=±{abs(xm):6.1f}")
