#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
子弹星团的真实双峰透镜模型
===========================

用文献给的实测参数，不是随手调的：

    主团   ρ0 = 2.6e-25 g/cm³,  rc = 210 kpc   （Markevitch+2004 King 剖面）
    子团   ρ0 = 2.0e-24 g/cm³,  rc = 54 kpc
    两峰间距  748 kpc（BCG 间距，Clowe+2006）

Clowe+2006 的探测显著性：主峰 12σ、子峰 8σ —— 两个峰都是硬的。

问题：
    标准做法是「两个势阱线性叠加」。它给出的两峰之间是什么形状？
    如果《构造》主张「两个凹陷之间还有膜张力耦合」，那多出来的东西
    会落在哪里、长什么样？

这个脚本只算**标准叠加**（作为基线）。基线算出来，才知道观测上
该去找什么形状的偏差。

用法：
    python bullet_lens_model.py
    python bullet_lens_model.py --plot ../figures/22_bullet_lens.png
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

KPC_IN_CM = 3.086e21

# Markevitch+2004 King 剖面参数
MAIN = dict(name="主团", rho0=2.6e-25, rc=210.0, x=-374.0)
SUB = dict(name="子团", rho0=2.0e-24, rc=54.0, x=+374.0)

SIGMA_CRIT = 0.70   # g/cm²，z=0.296 透镜、源在 z~1 量级


def sigma_surface(R, rho0, rc):
    """King 剖面的投影面密度 (g/cm²)。R, rc 单位 kpc。"""
    rc_cm = rc * KPC_IN_CM
    R_cm = R * KPC_IN_CM
    return 2.0 * rho0 * rc_cm / math.sqrt(1.0 + (R_cm / rc_cm) ** 2)


def kappa(R, lens):
    return sigma_surface(R, lens["rho0"], lens["rc"]) / SIGMA_CRIT


def enclosed_kappa(R_max, lens, n=4000):
    """∫0^R κ(R') 2πR' dR' / (πR²) —— 用于算偏折和剪切。"""
    if R_max <= 0:
        return 0.0
    dR = R_max / n
    tot = 0.0
    for i in range(n):
        Ri = (i + 0.5) * dR
        tot += kappa(Ri, lens) * Ri * dR
    return 2.0 * tot / (R_max * R_max)


def axially_symmetric(R, lens):
    """轴对称透镜的 κ 和 γ。γ = <κ> - κ。"""
    k = kappa(R, lens)
    kbar = enclosed_kappa(R, lens)
    return k, max(kbar - k, 0.0)


def two_lens_profile(x):
    """沿连线（y=0）的总 κ 和 |γ|（标量叠加，示意用）。"""
    k_tot = g_tot = 0.0
    for L in (MAIN, SUB):
        R = abs(x - L["x"])
        k, g = axially_symmetric(R, L)
        k_tot += k
        g_tot += g
    return k_tot, g_tot


def build_report():
    out = []
    add = out.append
    add("=" * 78)
    add("子弹星团双峰透镜模型（文献参数，非手调）")
    add("=" * 78)
    add("主团 ρ0=%.1e g/cm³  rc=%.0f kpc   位置 x=%.0f kpc"
        % (MAIN["rho0"], MAIN["rc"], MAIN["x"]))
    add("子团 ρ0=%.1e g/cm³  rc=%.0f kpc   位置 x=%.0f kpc"
        % (SUB["rho0"], SUB["rc"], SUB["x"]))
    add("两峰间距 %.0f kpc    Σ_crit = %.2f g/cm²" % (SUB["x"] - MAIN["x"], SIGMA_CRIT))
    add("")

    # 核心 κ
    for L in (MAIN, SUB):
        add("  %s 核心 κ = %.2f" % (L["name"], kappa(0.0, L)))
    add("  （κ < 1 → 弱透镜为主；子团核心接近 1，对应观测到的强透镜弧）")
    add("")

    add("-" * 78)
    add("一、沿两峰连线的剖面")
    add("-" * 78)
    add("   x(kpc) |   κ(主)  |   κ(子)  |   κ 总  |   γ 总  | 位置")
    add("  --------+---------+---------+---------+---------+------")
    for x in (-500, -450, -400, -374, -300, -200, -100, 0,
              100, 200, 300, 374, 400, 450, 500):
        k1, _ = axially_symmetric(abs(x - MAIN["x"]), MAIN)
        k2, _ = axially_symmetric(abs(x - SUB["x"]), SUB)
        k_tot, g_tot = two_lens_profile(float(x))
        tag = ""
        if abs(x - MAIN["x"]) < 1:
            tag = "主峰"
        elif abs(x - SUB["x"]) < 1:
            tag = "子峰"
        elif abs(x) < 1:
            tag = "两峰中间"
        add("  %7d | %7.3f | %7.3f | %7.3f | %7.4f | %s"
            % (x, k1, k2, k_tot, g_tot, tag))

    add("")
    add("-" * 78)
    add("二、两峰之间（x = 0 附近）的形状")
    add("-" * 78)
    add("  标准模型（线性叠加）在这里给出的是一条**单调的鞍部**：")
    add("  从左峰下来、经过中点、再上右峰。没有局部凸起。")
    add("")
    add("  量一下这个鞍部有多平：")
    k_mid, g_mid = two_lens_profile(0.0)
    k_left, _ = two_lens_profile(-200.0)
    k_right, _ = two_lens_profile(200.0)
    add("    x=-200：κ=%.3f   x=0：κ=%.3f   x=+200：κ=%.3f" % (k_left, k_mid, k_right))
    add("    连线的 κ 在中点附近 ≈ %.3f，占总范围的比例不低——" % k_mid)
    add("    因为**主团的核心半径 210 kpc 比两峰间距的一半（374 kpc）还大**。")
    add("    也就是说：主团的势阱一直铺到子峰脚底下，两峰之间的区域")
    add("    本来就是主团的「内部」，不是真正的空地。")
    add("")
    add("  ★ 这一点很关键：所谓「两峰之间有鞍部」，鞍部并不深。")
    add("    主团的 rc=210 kpc 决定了它自己就填满了中间那块。")
    add("")

    add("-" * 78)
    add("三、《构造》的额外项该落在哪里")
    add("-" * 78)
    add("  如果两根绳缆滑脱后还留着张力耦合，那耦合只能作用在")
    add("  **两个凹陷都够强的地方** —— 也就是两峰之间那一段。")
    add("")
    add("  标准模型在这里给 κ ≈ %.2f、γ ≈ %.4f。" % (k_mid, g_mid))
    add("  额外的「脊」会表现为：")
    add("    (a) κ 在中点附近高于两条 King 剖面的和")
    add("    (b) 剪切方向在中点附近出现**额外扭转**（因为脊不是轴对称的）")
    add("")
    add("  (b) 比 (a) 好找：κ 的绝对值受 Σ_crit 和源红移影响，")
    add("  而剪切的方向是纯几何量，不依赖归一化。")
    add("")
    add("  ★ 所以检验方法：**量两峰之间那段的剪切方向**。")
    add("    标准模型预言它单调地从「沿连线」转到「垂直连线」再转回来；")
    add("    《构造》如果多一条脊，就会在这个单调过渡里插进一个**额外的折点**。")
    return "\n".join(out)


def make_plot(path):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    # 网格
    N = 500
    L = 900.0
    xs = np.linspace(-L, L, N)
    ys = np.linspace(-L, L, N)
    X, Y = np.meshgrid(xs, ys)

    def prof(R, lens):
        rc_cm = lens["rc"] * KPC_IN_CM
        return 2.0 * lens["rho0"] * rc_cm / np.sqrt(
            1.0 + (R * KPC_IN_CM / rc_cm) ** 2) / SIGMA_CRIT

    K = np.zeros_like(X)
    for lens in (MAIN, SUB):
        R = np.hypot(X - lens["x"], Y)
        K += prof(R, lens)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.4))

    # 左：二维 κ
    ax = axes[0]
    im = ax.imshow(K, extent=[-L, L, -L, L], cmap="inferno", origin="lower",
                   vmin=0, vmax=min(K.max(), 1.0))
    plt.colorbar(im, ax=ax, fraction=0.046, label="会聚 κ")
    ax.contour(X, Y, K, levels=[0.1, 0.2, 0.3, 0.4, 0.5],
               colors="cyan", linewidths=0.9, alpha=0.8)
    ax.plot([MAIN["x"], SUB["x"]], [0, 0], "w+", ms=14, mew=2.4)
    ax.annotate("主团 12σ", (MAIN["x"], 0), xytext=(MAIN["x"] - 40, 120),
                color="w", fontsize=9, ha="center",
                arrowprops=dict(arrowstyle="->", color="w", lw=1))
    ax.annotate("子团 8σ", (SUB["x"], 0), xytext=(SUB["x"] + 40, 120),
                color="w", fontsize=9, ha="center",
                arrowprops=dict(arrowstyle="->", color="w", lw=1))
    ax.set_xlabel("x (kpc)")
    ax.set_ylabel("y (kpc)")
    ax.set_title("子弹星团：两个实测质量峰的会聚图\n（King 剖面参数取自 Markevitch+2004，峰位取自 Clowe+2006）",
                 fontsize=10)

    # 右：沿连线剖面
    ax = axes[1]
    xline = np.linspace(-700, 700, 700)
    k_tot = np.array([two_lens_profile(float(x))[0] for x in xline])
    g_tot = np.array([two_lens_profile(float(x))[1] for x in xline])

    ax.plot(xline, k_tot, lw=2.4, color="#c53030", label="κ 总（线性叠加）")
    ax.plot(xline, g_tot, lw=2.0, color="#2b6cb0", ls="--", label="γ 总")
    ax.axvline(MAIN["x"], color="gray", ls=":", lw=1.2)
    ax.axvline(SUB["x"], color="gray", ls=":", lw=1.2)
    ax.axvline(0, color="#38a169", ls="-.", lw=1.4, alpha=0.8)
    ax.text(MAIN["x"], max(k_tot) * 0.95, " 主峰", fontsize=9, color="gray")
    ax.text(SUB["x"], max(k_tot) * 0.95, " 子峰", fontsize=9, color="gray")
    ax.text(6, max(k_tot) * 0.72, "中点", fontsize=9, color="#38a169")

    # 标出鞍部
    k_mid = two_lens_profile(0.0)[0]
    ax.axhline(k_mid, color="#38a169", ls=":", lw=1.0, alpha=0.6)
    ax.annotate("鞍部 κ≈%.2f" % k_mid, xy=(0, k_mid), xytext=(120, k_mid - 0.12),
                fontsize=9, color="#38a169",
                arrowprops=dict(arrowstyle="->", color="#38a169", lw=1))

    ax.set_xlabel("沿两峰连线的位置 x (kpc)")
    ax.set_ylabel("κ / γ")
    ax.set_title("沿连线剖面：两峰之间是一条浅鞍部\n（主团 rc=210 kpc 自己就填满了中间）",
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)

    fig.suptitle("子弹星团的真实双峰透镜 —— 基线模型", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="子弹星团双峰透镜模型")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    text = build_report()
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
        print("report written: %s" % args.out)
    else:
        print(text)

    if args.plot:
        make_plot(args.plot)


if __name__ == "__main__":
    main()
