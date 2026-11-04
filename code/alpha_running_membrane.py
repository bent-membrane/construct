# -*- coding: utf-8 -*-
"""
α 随尺度跑动 → 膜参数随尺度变化

关系（下篇 4.9 解出来的）：
    d    = l_P · sqrt(4·pi·α)
    T    = (F_P/l_P) / sqrt(4·pi·α)
    T/d  = p_P / (4·pi·α)

α 的跑动：两端点插值
    α⁻¹(0)    = 137.035999
    α⁻¹(M_Z)  = 127.952
形状用 ln(1+(Q/Q0)²)，Q0=0.1 GeV（强子阈值附近）。
**这是示意模型，不是精确 QED 单圈**——只保证两端点对、形状单调。
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --plot <path> 才保存图片（符合仓库约定：直接运行只打印数值）
PLOT = None
if "--plot" in sys.argv:
    i = sys.argv.index("--plot")
    if i + 1 < len(sys.argv):
        PLOT = sys.argv[i + 1]

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

# ---------- 常数 ----------
c = 299792458.0
G = 6.67430e-11
hbar = 1.054571817e-34
lP = np.sqrt(hbar * G / c**3)      # 1.616255e-35 m
FP = c**4 / G                       # 1.2103e44 N
pP = c**7 / (hbar * G**2)           # 4.633e113 Pa

# ---------- α 跑动（两端点插值） ----------
A0_INV = 137.035999
AZ_INV = 127.952
MZ = 91.1876       # GeV
Q0 = 0.1           # GeV

def f_run(Q):
    return np.log1p((Q / Q0) ** 2) / np.log1p((MZ / Q0) ** 2)

def alpha_inv(Q):
    return A0_INV - (A0_INV - AZ_INV) * f_run(Q)

def alpha(Q):
    return 1.0 / alpha_inv(Q)

# ---------- 扫描 ----------
Qs = np.logspace(-3, 3, 400)        # 1 MeV → 1 TeV
a = alpha(Qs)

d = lP * np.sqrt(4 * np.pi * a)
T = (FP / lP) / np.sqrt(4 * np.pi * a)
P = pP / (4 * np.pi * a)

d0, T0, P0 = d[0], T[0], P[0]

# ---------- 关键点 ----------
i0 = 0
iZ = int(np.argmin(np.abs(Qs - MZ)))

print("=" * 62)
print("α 跑动 → 膜参数")
print("=" * 62)
print(f"{'Q':>12} {'alpha^-1':>12} {'d/lP':>10} {'T (N/m)':>12} {'T/d (Pa)':>12}")
for i in [0, iZ, len(Qs) - 1]:
    print(f"{Qs[i]:10.3e} GeV {alpha_inv(Qs[i]):12.3f} {d[i]/lP:10.5f} "
          f"{T[i]:12.4e} {P[i]:12.4e}")
print()
print(f"d/lP        : {d0/lP:.5f}  ->  {d[iZ]/lP:.5f}   ({(d[iZ]/d0-1)*100:+.2f}%)")
print(f"T           : {T0:.4e} ->  {T[iZ]:.4e}  ({(T[iZ]/T0-1)*100:+.2f}%)")
print(f"T/d         : {P0:.4e} ->  {P[iZ]:.4e}  ({(P[iZ]/P0-1)*100:+.2f}%)")
print()
print("方向：α 随能量变大 → d 变厚、T 变小、T/d 变小（膜变软）")

# ---------- 画图 ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

ax = axes[0]
ax.semilogx(Qs, alpha_inv(Qs), lw=2, color="#C0392B")
ax.axvline(MZ, ls="--", lw=1, color="gray")
ax.axhline(AZ_INV, ls=":", lw=1, color="gray")
ax.set_xlabel("能量尺度 Q (GeV)")
ax.set_ylabel(r"$\alpha^{-1}$")
ax.set_title(r"① $\alpha$ 随尺度跑动" + "\n" + r"(137.036 → 127.95)")
ax.grid(alpha=0.3)

ax = axes[1]
ax.semilogx(Qs, d / lP, lw=2, color="#2980B9")
ax.axvline(MZ, ls="--", lw=1, color="gray")
ax.axhline(d0 / lP, ls=":", lw=1, color="gray")
ax.set_xlabel("能量尺度 Q (GeV)")
ax.set_ylabel(r"$d\,/\,l_P$")
ax.set_title(r"② 膜厚 = $l_P\sqrt{4\pi\alpha}$" + "\n" + f"({d0/lP:.4f} → {d[iZ]/lP:.4f} lP)")
ax.grid(alpha=0.3)

ax = axes[2]
ax.semilogx(Qs, (T / T0 - 1) * 100, lw=2, label=r"$T$（张力）", color="#27AE60")
ax.semilogx(Qs, (P / P0 - 1) * 100, lw=2, label=r"$T/d$（内应力）", color="#8E44AD")
ax.axvline(MZ, ls="--", lw=1, color="gray")
ax.axhline(0, ls=":", lw=1, color="gray")
ax.set_xlabel("能量尺度 Q (GeV)")
ax.set_ylabel("相对 Q=0 的变化 (%)")
ax.set_title("③ 张力与内应力随尺度变化")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

fig.suptitle(r"$\alpha$ 跑动 → $d=l_P\sqrt{4\pi\alpha}$：能量越高，膜越厚、越软", fontsize=12)
fig.tight_layout()
if PLOT:
    fig.savefig(PLOT, dpi=160)
    print("\nSaved:", PLOT)
else:
    print("\n(加 --plot <path> 保存图片)")
