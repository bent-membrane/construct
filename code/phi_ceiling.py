#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
φ 有没有天花板？—— 《构造》核/雾相位轨迹的判据
=================================================

⚠ 降级 + 已判死（2026-10-29）：
   ① 骨架来自 four_clusters.py，时间/速度只作相位刻度；
   ② φ 无天花板——400 组随机扫描 φ_max 全部 > 87°，
      「任何并合系统都不该超过 X°」这个形式的判据不成立；
   ③ 残余判据（φ>50° ⇒ 核在雾后 + ΔS 小）不可执行：「核在雾后」需独立相位判定，
      单系统测不了；四个系统最大 41°，判据区为空；
   ④ U90 这条线关闭。

模型与 four_clusters.py 完全配套（同一套骨架、同一套参数）：

    核：阻尼振荡      r(t) = (v0/wd)·exp(-γt)·sin(wd·t)
    雾：只进不退      s(t) = 震落后减速的位移
    纵向：            dSx = r − s
    横向：回头前 b 线性长；回头后按**回程距离**幂律长
                      dSy = b + k·(回程距离)^α

先说结论方向（跑完验证）：

    φ 的峰值出现在 **核穿过雾的那一瞬**（dSx → 0）。
    那一瞬 tan φ = dSy/dSx → ∞，所以 φ 能到 90°，**不是硬上限**。

    但那个窗口很窄。核越过雾的相对速度 ~2400 km/s，
    要 |dSx| < 100 kpc，窗口只有 ~80 Myr。

    所以正确的说法不是「φ 超不过 50°」，而是：
        **φ > 50° 的系统必然同时满足：核在雾后面 + ΔS 很小**
    这是一个可以被单个反例杀死的判据。

本脚本做三件事：
    1. 用与 18 图同一套参数，看 φ 怎么冲顶又回落
    2. 扫参数，看 φ_max 的分布（是不是真的能到 90）
    3. 给出可证伪的判据：φ 与 ΔS / 核方位的联合约束

用法：
    python phi_ceiling.py
    python phi_ceiling.py --plot ../figures/19_phi_ceiling.png
"""

import argparse
import math
import random
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

KMPS_MYR_TO_KPC = 1.0227e-3


# ---------------------------------------------------------------- 与 18 图同一套骨架
def core_r(t, v0, wd, gamma):
    if t <= 0:
        return 0.0
    return (v0 / wd) * math.exp(-gamma * t) * math.sin(wd * t)


def fog_s(t, v0, vinf, tau, n=300):
    if t <= 0:
        return 0.0
    dt = t / n
    total = 0.0
    for i in range(n):
        tt = (i + 0.5) * dt
        total += (vinf + (v0 - vinf) * math.exp(-tt / tau)) * dt
    return total


def traj(t, p):
    """返回 dSx, dSy, ΔS, φ观测(锐角), φ带符号。"""
    r = core_r(t, p["v0"], p["wd"], p["gamma"])
    s = fog_s(t, p["v0"], p["v0"] / p["q"], p["tau"])
    dSx = (r - s) * KMPS_MYR_TO_KPC

    t_turn = math.pi / (2 * p["wd"])
    if t <= t_turn:
        dSy = p["b"] * (t / t_turn)
    else:
        # 甩动按「回程距离」起步 —— 与 four_clusters.py 一致，避免接缝
        r_kpc = r * KMPS_MYR_TO_KPC
        r_max = core_r(t_turn, p["v0"], p["wd"], p["gamma"]) * KMPS_MYR_TO_KPC
        back = max(0.0, r_max - r_kpc)
        dSy = p["b"] + p["k_swing"] * (back ** p["alpha"])

    dS = math.hypot(dSx, dSy)
    phi_s = math.degrees(math.atan2(dSy, dSx)) if abs(dSx) > 1e-9 else 90.0
    phi = phi_s if phi_s <= 90 else 180 - phi_s
    return dSx, dSy, dS, phi, phi_s


def default_params():
    """与 four_clusters.default_params() 对齐。"""
    return dict(
        v0=2650.0,
        q=10.0,
        tau=40.0,
        b=30.5,
        k_swing=3.2e-4,
        alpha=2.28,
        wd=1.0 / 220.0,
        gamma=1.0 / 600.0,
    )


def scan_curve(p, t_max=2000.0, n=1200):
    """与 18 图同一窗口（2000 Myr）。超过一个半周期后幂律会外推失控：
    回程距离继续变大会让 dSy 涨到非物理量级，曲线在图上冲出象限。"""
    return [(t_max * i / n, traj(t_max * i / n, p)) for i in range(1, n + 1)]


def ceiling(p):
    """返回 (φ_max, t_at_max, 高φ时间占比, 周期)"""
    curve = scan_curve(p)
    phi_max, t_max = -1.0, 0.0
    for t, (_, _, _, phi, _) in curve:
        if phi > phi_max:
            phi_max, t_max = phi, t
    n_hi = sum(1 for _, (_, _, _, phi, _) in curve if phi > 50.0)
    return phi_max, t_max, n_hi / len(curve), 2 * math.pi / p["wd"]


def build_report(p):
    out = []
    add = out.append
    add("=" * 74)
    add("φ 天花板检测 — 《构造》核/雾轨迹（与 18 图同一套参数）")
    add("=" * 74)
    add("v0=%.0f  q=%.1f  tau=%.0f  b=%.1f  k_swing=%.2e  alpha=%.2f"
        % (p["v0"], p["q"], p["tau"], p["b"], p["k_swing"], p["alpha"]))
    add("周期 = %.0f Myr   阻尼时间 = %.0f Myr"
        % (2 * math.pi / p["wd"], 1.0 / p["gamma"]))
    add("")

    curve = scan_curve(p)
    zc = []
    prev = None
    for t, (dSx, _, _, _, _) in curve:
        if prev is not None and prev[1] * dSx < 0:
            zc.append(prev[0])
        prev = (t, dSx)
    add("dSx 过零（核穿过雾）时刻：" +
        (", ".join("%.0f Myr" % z for z in zc) if zc else "无"))
    add("")

    add("-" * 74)
    add("穿过雾前后 300 Myr：φ 怎么冲上去又回落")
    add("-" * 74)
    add("   t(Myr) |  dSx   |  dSy   |   ΔS   |  φ(°)  | 核位置")
    add("  ---------+--------+--------+--------+--------+--------")
    t_ref = zc[0] if zc else 500.0
    for dt in (-300, -200, -120, -60, -20, 0, 20, 60, 120, 200, 300):
        t = max(t_ref + dt, 1.0)
        dSx, dSy, dS, phi, _ = traj(t, p)
        where = "雾前" if dSx > 0 else "雾后"
        add("  %8.0f | %6.0f | %6.0f | %6.0f | %6.1f | %s"
            % (t, dSx, dSy, dS, phi, where))

    add("")
    add("-" * 74)
    add("参数扫描：φ_max 能到多少（其余参数固定为 18 图那套）")
    add("-" * 74)
    add("  参数        取值         φ_max     高φ窗口占比")
    add("  -----------+-----------+----------+------------")
    grid = {
        "q": [3, 5, 10, 20, 40],
        "tau": [20, 40, 80, 150],
        "k_swing": [1e-4, 3.2e-4, 1e-3, 3e-3, 1e-2],
        "alpha": [1.5, 2.0, 2.28, 3.0, 4.0],
        "period": [800, 1100, 1400, 2000, 2800],
        "damp": [300, 600, 1200, 3000],
    }
    for key, vals in grid.items():
        for v in vals:
            q = dict(p)
            if key == "period":
                q["wd"] = 2 * math.pi / v
            elif key == "damp":
                q["gamma"] = 1.0 / v
            else:
                q[key] = v
            pm, tm, frac, period = ceiling(q)
            add("  %-10s %11s  %6.1f°    %6.1f%%" % (key, v, pm, 100 * frac))
    add("")

    add("-" * 74)
    add("随机扫 400 组参数：φ_max 的分布")
    add("-" * 74)
    random.seed(7)
    pms = []
    for _ in range(400):
        q = dict(p)
        q["q"] = random.uniform(2, 50)
        q["tau"] = random.uniform(10, 200)
        q["k_swing"] = 10 ** random.uniform(-4.5, -2)
        q["alpha"] = random.uniform(1.2, 4.0)
        q["wd"] = 2 * math.pi / random.uniform(700, 3000)
        q["gamma"] = 1.0 / random.uniform(200, 4000)
        pm, tm, frac, _ = ceiling(q)
        pms.append(pm)
    pms.sort()
    n = len(pms)
    add("  最小 %.1f°   10%% %.1f°   中位 %.1f°   90%% %.1f°   最大 %.1f°"
        % (pms[0], pms[n // 10], pms[n // 2], pms[9 * n // 10], pms[-1]))
    add("  φ_max > 50° 的比例：%.0f%%" % (100 * sum(1 for x in pms if x > 50) / n))
    add("  φ_max > 70° 的比例：%.0f%%" % (100 * sum(1 for x in pms if x > 70) / n))
    add("")

    add("-" * 74)
    add("结论")
    add("-" * 74)
    add("  φ **没有硬上限**：核穿过雾的那一瞬 dSx→0，φ 数学上可以逼近 90°。")
    add("  但高 φ 的窗口很窄（核相对雾 ~2400 km/s，|dSx|<100 kpc 只占 ~80 Myr），")
    add("  在 1-2 Gyr 的周期里是小概率事件。")
    add("")
    add("  可证伪的判据因此不是「φ 有上限」，而是：")
    add("    ★ 凡是 φ > 50° 的系统，必须同时满足")
    add("        (a) 核在雾的后面（反向分离）")
    add("        (b) ΔS 不大（|dSx| 小，通常 < 250 kpc）")
    add("    如果找到一个 φ > 50° 但 ΔS > 400 kpc、且核在雾前面的系统，")
    add("    这条轨迹就错了。")
    return "\n".join(out)


def make_plot(p, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.2))

    # 左：φ(t) 曲线族（按甩动强度区分）
    ax = axes[0]
    base = scan_curve(p, n=900)
    ts = [t for t, _ in base]
    ax.plot(ts, [v[3] for _, v in base], lw=2.6, color="#c53030",
            label="基准（k=%.1e）" % p["k_swing"])

    for k, color in ((3e-5, "#a0aec0"), (3e-4, "#718096"), (3e-3, "#2d3748")):
        q = dict(p); q["k_swing"] = k
        c = scan_curve(q, n=900)
        ax.plot([t for t, _ in c], [v[3] for _, v in c], lw=1.4, ls="--",
                color=color, label="k = %.0e" % k)

    ax.axhspan(50, 90, color="#dd6b20", alpha=0.12, label="φ > 50° 窗口")
    ax.axhline(41, color="#2b6cb0", ls=":", lw=1.4, label="观测最高 41°（ZwCl）")
    ax.axhline(90, color="black", ls="-.", lw=1.0, alpha=0.5, label="理论可达 90°")
    ax.set_xlabel("碰撞后时间 t (Myr)")
    ax.set_ylabel("分离角 φ (°)")
    ax.set_title("φ 能冲多高：核穿过雾时逼近 90°，但窗口很窄")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.25)
    ax.set_ylim(0, 95)

    # 右：φ vs ΔS 联合约束
    ax = axes[1]
    for k, color, alpha in ((3e-5, "#a0aec0", 0.55), (3e-4, "#718096", 0.7),
                            (3.2e-4, "#c53030", 0.9), (3e-3, "#2d3748", 0.75)):
        q = dict(p); q["k_swing"] = k
        c = scan_curve(q, n=900)
        ax.scatter([v[2] for _, v in c], [v[3] for _, v in c],
                   s=9, alpha=alpha, color=color, label="k = %.0e" % k)

    ax.axhspan(50, 95, color="#dd6b20", alpha=0.10)
    ax.text(360, 62, "判据区：\nφ>50° 必须 ΔS 小", fontsize=9, color="#dd6b20",
            ha="center", fontweight="bold")

    # 观测框（不是中心点）：模型穿的是框的边角，画中心会假装匹配
    obs = [("MACS J0025", 170, 370, 2.0, 5.0, "#718096"),
           ("Bullet", 100, 162, 11.8, 14.1, "#c53030"),
           ("A520", 368, 630, 21.5, 39.0, "#dd6b20"),
           ("ZwCl 0008", 173, 312, 32.0, 41.0, "#2b6cb0")]
    for name, x0, x1, y0, y1, color in obs:
        ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=True,
                                   alpha=0.20, color=color, ec=color, lw=1.4,
                                   zorder=3))
        ax.annotate(name, xy=(0.5 * (x0 + x1), y1), xytext=(0, 4),
                    textcoords="offset points", ha="center", fontsize=8.5,
                    color=color, fontweight="bold")

    ax.set_xlabel("核—雾分离距离 ΔS (kpc)")
    ax.set_ylabel("分离角 φ (°)")
    ax.set_title("联合约束：φ 与 ΔS 不能同时大")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.25)
    ax.set_xlim(0, 700)
    ax.set_ylim(0, 95)

    fig.suptitle("《构造》φ 天花板判据 —— 不是硬上限，是窄窗口", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="φ 天花板判据（与 18 图配套）")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    p = default_params()
    text = build_report(p)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
        print("report written: %s" % args.out)
    else:
        print(text)

    if args.plot:
        make_plot(p, args.plot)


if __name__ == "__main__":
    main()
