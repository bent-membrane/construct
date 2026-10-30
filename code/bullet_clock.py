# -*- coding: utf-8 -*-
"""
图 04：给双曲线配一块表 —— 子弹星团
X 射线激波 / 碰撞时间 当"时钟"，看能不能把速度钉死。
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

OUT = os.path.dirname(os.path.abspath(__file__))

YR = 3.1558e7          # s
KPC = 3.0857e19        # m
MPC = 1e3 * KPC
KMS = 1e3              # m/s

# 子弹星团（1E 0657-56）
SEP_KPC = 720.0                    # 两个子团在天空上的分离（投影）
V_LO, V_HI = 2700.0, 4700.0        # 文献里激波速度的取值范围 km/s

sep = SEP_KPC * KPC
v_lo, v_hi = V_LO * KMS, V_HI * KMS

tau_max = sep / v_lo / YR          # 用慢端算 => 最长可能时间
tau_min = sep / v_hi / YR          # 用快端算 => 最短可能时间
print(f"分离 {SEP_KPC:.0f} kpc")
print(f"激波速度 {V_LO:.0f}-{V_HI:.0f} km/s")
print(f"-> 碰撞至今 {tau_min/1e8:.2f}-{tau_max/1e8:.2f} 亿年")

TAU_REF = 2.35e8                   # 取中间值做示范

fig, axes = plt.subplots(1, 2, figsize=(14, 5.4))

# ---------------------------------------------------------------- A：给时间
ax = axes[0]
R = np.linspace(0.15, 2.2, 300) * MPC
for tau_yr, c, ls in [
    (1.5e8, "#8e44ad", ":"),
    (TAU_REF, "#c0392b", "-"),
    (3.0e8, "#2980b9", "--"),
]:
    v = R / (tau_yr * YR) / KMS
    ax.plot(R / MPC, v, color=c, ls=ls, lw=2, label=f"如果年龄是 {tau_yr/1e8:.2f} 亿年")

ax.axhspan(V_LO, V_HI, color="#f39c12", alpha=0.25,
           label=f"X 射线激波测到 {V_LO:.0f}–{V_HI:.0f} km/s")

# 自洽区
r_lo = v_lo * TAU_REF * YR / MPC
r_hi = v_hi * TAU_REF * YR / MPC
ax.axvspan(r_lo, r_hi, color="#27ae60", alpha=0.12)
ax.plot([r_lo, r_hi], [V_LO, V_HI], "o", color="#27ae60", ms=7, zorder=5)
ax.annotate(f"半径 {r_lo:.2f}–{r_hi:.2f} Mpc\n两条独立测量在这里碰头",
            xy=((r_lo + r_hi) / 2, (V_LO + V_HI) / 2),
            xytext=(0.55, 7600), fontsize=10,
            arrowprops=dict(arrowstyle="->", color="#27ae60"))

ax.set_xlabel("观测到的雾化层半径（Mpc）")
ax.set_ylabel("喷射速度（km/s）")
ax.set_title("A. 年龄当尺子：速度 = 半径 ÷ 年龄", fontsize=12)
ax.set_xlim(0.15, 2.2)
ax.set_ylim(0, 10000)
ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.grid(alpha=0.2)

# ---------------------------------------------------------------- B：给速度
ax2 = axes[1]
tau = np.logspace(7.4, 9.2, 300)      # yr
for R_Mpc, c, ls in [
    (0.5, "#8e44ad", ":"),
    (1.0, "#c0392b", "-"),
    (1.5, "#2980b9", "--"),
]:
    v = (R_Mpc * MPC) / (tau * YR) / KMS
    ax2.plot(tau, v, color=c, ls=ls, lw=2, label=f"如果半径是 {R_Mpc:.1f} Mpc")

ax2.axhspan(V_LO, V_HI, color="#f39c12", alpha=0.25, label="激波测到的速度")
ax2.axvspan(tau_min, tau_max, color="#27ae60", alpha=0.12, label="碰撞时间（由分离÷速度推出）")
ax2.plot([TAU_REF], [1.0 * MPC / (TAU_REF * YR) / KMS], "o",
         ms=10, mfc="none", mec="#27ae60", mew=2.5, zorder=6)

ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("喷流年龄（年）")
ax2.set_ylabel("喷射速度（km/s）")
ax2.set_title("B. 原来的双曲线：加一块表，就落成一个点", fontsize=12)
ax2.legend(frameon=False, fontsize=8.5, loc="lower left")
ax2.grid(alpha=0.2, which="both")

plt.tight_layout()
p = os.path.join(OUT, "04_bullet_clock.png")
plt.savefig(p, dpi=170, bbox_inches="tight")
plt.close()
print("saved", p)

# ---------------------------------------------------------------- 数
v_ref = 1.0 * MPC / (TAU_REF * YR) / KMS
print(f"\n半径 1.0 Mpc + 年龄 {TAU_REF/1e8:.2f} 亿年 -> 速度 {v_ref:.0f} km/s")
print(f"激波给的独立范围          -> {V_LO:.0f}-{V_HI:.0f} km/s")
print(f"落在里面吗: {V_LO <= v_ref <= V_HI}")
print(f"\n反过来说：若速度取 3000 km/s，1 Mpc 对应年龄 {1.0*MPC/(3000*KMS)/YR/1e8:.2f} 亿年")
