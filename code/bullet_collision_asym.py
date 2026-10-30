# -*- coding: utf-8 -*-
"""
图 06：子弹星团（不对称版）
主团 : 子弹团 = 10 : 1
  - 主团：大、慢、气体几乎不被撞停
  - 子弹团：小、快、气体被撞停，雾照穿
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

KPC_KM = 3.0857e16
YR_S = 3.1558e7

# 主团（左，慢） / 子弹团（右，快）
VA = 470.0          # km/s，向右
VB = -3377.0        # km/s，向左
TAU_A = 200e6       # 主团气体几乎不减速
TAU_B = 150e6       # 子弹团气体被拖慢
XA0 = -24.0         # kpc
XB0 = 240.0         # kpc
T_END = 250e6

T_C = (XB0 - XA0) * KPC_KM / (VA - VB) / YR_S
print(f"碰撞时刻 t_c = {T_C/1e6:.1f} Myr")
print(f"相对速度 {VA-VB:.0f} km/s")


def pos(x0, v, tau, t):
    if t <= T_C:
        return x0 + v * t * YR_S / KPC_KM
    xc = x0 + v * T_C * YR_S / KPC_KM
    dx = v * tau * YR_S / KPC_KM * (1 - np.exp(-(t - T_C) / tau))
    return xc + dx


def gas(t):
    return pos(XA0, VA, TAU_A, t), pos(XB0, VB, TAU_B, t)


def mist(t):
    """雾：无碰撞，始终匀速（不受碰撞影响）。"""
    return (XA0 + VA * t * YR_S / KPC_KM,
            XB0 + VB * t * YR_S / KPC_KM)


# ---------------------------------------------------------------- 画布
fig, ax = plt.subplots(figsize=(10, 5.2), facecolor="#05060a")
ax.set_facecolor("#05060a")
ax.set_xlim(-800, 500)
ax.set_ylim(-190, 190)
ax.set_aspect("equal")
ax.set_xlabel("膜上位置（kpc）", color="#bbbbbb")
ax.set_title("子弹星团（不对称）：大团慢、小团快，气体撞停、雾穿透",
             color="white", fontsize=12.5, pad=10)
ax.tick_params(colors="#8899aa")
for s in ax.spines.values():
    s.set_color("#223344")

rng = np.random.default_rng(5)
# 主团：大而弥散；子弹团：小而密
NA, NB = 3200, 900
dA_x, dA_y = rng.normal(0, 55, NA), rng.normal(0, 55, NA)
dB_x, dB_y = rng.normal(0, 22, NB), rng.normal(0, 22, NB)

mist_A = ax.scatter([], [], s=2.4, c="#3fa9f5", alpha=0.14, linewidths=0)
mist_B = ax.scatter([], [], s=2.6, c="#5ec8ff", alpha=0.16, linewidths=0)

gA = ax.scatter([], [], s=300, c="#ff7b00", alpha=0.16, linewidths=0)
gB = ax.scatter([], [], s=200, c="#ff7b00", alpha=0.20, linewidths=0)
cA = ax.scatter([], [], s=70, c="#ffd166", alpha=0.95, linewidths=0)
cB = ax.scatter([], [], s=55, c="#ffd166", alpha=0.95, linewidths=0)

txt_t = ax.text(0.02, 0.95, "", transform=ax.transAxes, color="white",
                fontsize=16, va="top", family="monospace")
txt_d = ax.text(0.02, 0.06, "", transform=ax.transAxes, color="#9fb3c8",
                fontsize=9.5, va="bottom")

ax.scatter([], [], s=45, c="#ffd166", label="气体（中心点）")
ax.scatter([], [], s=45, c="#3fa9f5", alpha=0.5, label="雾（投影）")
ax.legend(loc="upper right", frameon=False, labelcolor="white", fontsize=9)

FRAMES = np.arange(0, T_END + 1, 5e6)


def update(t):
    xm_a, xm_b = mist(t)
    xg_a, xg_b = gas(t)

    mist_A.set_offsets(np.column_stack([xm_a + dA_x, dA_y]))
    mist_B.set_offsets(np.column_stack([xm_b + dB_x, dB_y]))
    gA.set_offsets([[xg_a, 0.0]])
    gB.set_offsets([[xg_b, 0.0]])
    cA.set_offsets([[xg_a, 0.0]])
    cB.set_offsets([[xg_b, 0.0]])

    txt_t.set_text(f"t = {t/1e6:5.0f} Myr")
    txt_d.set_text(
        f"雾间距 {abs(xm_a-xm_b):5.0f} kpc    "
        f"气体间距 {abs(xg_a-xg_b):5.0f} kpc    "
        f"子弹团 气体与雾 分开 {abs(xg_b-xm_b):5.0f} kpc"
    )
    return mist_A, mist_B, gA, gB, cA, cB, txt_t, txt_d


ani = FuncAnimation(fig, update, frames=FRAMES, interval=120, blit=False)
p = os.path.join(OUT, "06_bullet_collision_asym.gif")
ani.save(p, writer=PillowWriter(fps=8), dpi=110)
plt.close()
print("saved", p, f"{os.path.getsize(p)/1e6:.2f} MB")

print("\n关键帧：")
for t in [0, 50e6, T_C, 150e6, 250e6]:
    xm_a, xm_b = mist(t)
    xg_a, xg_b = gas(t)
    print(f"  t={t/1e6:5.1f}  主团雾={xm_a:7.1f} 主团气={xg_a:7.1f} | "
          f"子弹雾={xm_b:7.1f} 子弹气={xg_b:7.1f}")
