#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
四星系统一建模 —— 《构造》核/雾相位轨迹
=========================================

⚠ 降级（2026-10-29）：「阻尼振荡 r(t)」「振荡周期」「阻尼时间」只作**相位刻度**，
   不作物理预测。依据 2026-10-26 定「按几何拓扑考虑，不按质量/时间/速度考虑」；
   四个系统落在同一条相线上，靠的是 φ–ΔS 的几何位置，周期 1.4 Gyr 只是拟合刻度。

目标：不是给每个星系单独调参，而是让四个并合系统落在
「一根绳子解开、滑出去、被膜张力拉回来」的**同一条相位轨迹**上。

四个系统的观测（10-25 文献取值，只标量级）：

    系统              φ (°)          ΔS (kpc)        记忆里的相位判断
    MACS J0025.4    2 – 5          170 – 370        刚穿
    Bullet 1E0657   11.8 / 14.1    100 – 162        刚穿末
    A520            21.5 – 39      368 – 630        中间
    ZwCl 0008.8     32 – 41        173 – 312        回程

模型的三条腿：

  1. 核（断点）不是一辈子往外飞。滑出去以后，膜张力会把它往回拉。
     所以在主团系里，核的位移是**阻尼振荡**，不是匀速：

         r(t) = (v0 / wd) * exp(-γt) * sin(wd * t)

  2. 雾（尘）被震落后留在原地，被主团介质拖慢，**只进不退**：

         s(t) = 震落后减速的位移（和 bullet_core_fog.py 的 fog_s 同一套）

  3. 核—雾分离拆成纵向和横向两截：

         dSx(t) = r(t) - s(t)         纵向（沿碰撞轴）
         dSy(t) = b * (1 - exp(-t/τy))  横向（震落时甩出来的，随时间长大）

         φ(t) = atan(dSy / dSx)

  于是 φ 大有两种合法来源：
    (a) 纵向被拉没了（核回头、雾还在往前）→ dSx 变小
    (b) 横向被甩出来了（回程摇摆）        → dSy 变大

  —— 这正是作者 10-25 拍板的那句：「回来就是摇摆，会震出去也可以吸回来」

判据：
  把四个系统的观测框画在 φ–ΔS 平面上，看模型曲线能不能同时穿过四个框。
  **穿不过的地方，就是模型要改的地方，不是要调参遮掉的地方。**

用法：
    python four_clusters.py
    python four_clusters.py --plot ../figures/18_four_clusters.png
    python four_clusters.py --sweep          # 扫参数，看轨迹族
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

KMPS_MYR_TO_KPC = 1.0227e-3

# ---------------------------------------------------------------- 观测
# 数值取自 10-25 文献核查；区间是量法/文献差异，不是误差棒
CLUSTERS = [
    dict(name="MACS J0025", phi_lo=2.0,  phi_hi=5.0,   dS_lo=170, dS_hi=370,
         phase="刚穿",      color="#718096"),
    dict(name="Bullet 1E0657", phi_lo=11.8, phi_hi=14.1, dS_lo=100, dS_hi=162,
         phase="刚穿末",    color="#c53030"),
    dict(name="A520", phi_lo=21.5, phi_hi=39.0, dS_lo=368, dS_hi=630,
         phase="中间",      color="#dd6b20"),
    dict(name="ZwCl 0008", phi_lo=32.0, phi_hi=41.0, dS_lo=173, dS_hi=312,
         phase="回程",      color="#2b6cb0"),
]


# ---------------------------------------------------------------- 模型
def core_r(t, v0, wd, gamma):
    """核在主团系里的位移（km/s·Myr）。阻尼振荡，先出去再被拉回。"""
    if t <= 0:
        return 0.0
    return (v0 / wd) * math.exp(-gamma * t) * math.sin(wd * t)


def core_speed(t, v0, wd, gamma):
    """核的瞬时速度。过零 = 开始回头。"""
    if t <= 0:
        return v0
    e = math.exp(-gamma * t)
    return v0 * e * (math.cos(wd * t) - (gamma / wd) * math.sin(wd * t))


def fog_s(t, v0, vinf, tau, model="exp", n=600):
    """雾的位移（km/s·Myr）。只进不退。"""
    if t <= 0:
        return 0.0
    dt = t / n
    total = 0.0
    for i in range(n):
        tt = (i + 0.5) * dt
        if model == "exp":
            v = vinf + (v0 - vinf) * math.exp(-tt / tau)
        elif model == "quad":
            v = vinf + (v0 - vinf) / (1.0 + tt / tau)
        else:  # const
            v = vinf if (tau <= 0 or tt >= tau) else v0 - (v0 - vinf) * tt / tau
        total += v * dt
    return total


def trajectory(t, p):
    """返回任意时刻的 dSx / dSy / φ / ΔS。"""
    r = core_r(t, p["v0"], p["wd"], p["gamma"])
    s = fog_s(t, p["v0"], p["v0"] / p["q"], p["tau"], p.get("brake", "exp"))
    dSx = (r - s) * KMPS_MYR_TO_KPC

    # 横向两段式：回头之前只是碰撞参数留的残余，回头之后才甩起来
    t_turn = math.pi / (2 * p["wd"])
    if t <= t_turn:
        dSy = p["b"] * (t / t_turn)
    elif p.get("swing_mode") == "power":
        # 甩动幅度不饱和，按「回程距离」幂律长。
        # 不能用 |dSx|：回头那一刻 |dSx| 已有 ~190 kpc，会从 b 直接跳到
        # 190^alpha，图上出现一条竖线（这是接缝，不是物理）。
        r_kpc = r * KMPS_MYR_TO_KPC
        r_max = core_r(t_turn, p["v0"], p["wd"], p["gamma"]) * KMPS_MYR_TO_KPC
        back = max(0.0, r_max - r_kpc)
        dSy = p["b"] + p["k_swing"] * (back ** p["alpha"])
    else:
        dSy = p["b"] + p.get("b_swing", 0.0) * (
            1.0 - math.exp(-(t - t_turn) / p["tau_y"]))

    dS = math.hypot(dSx, dSy)
    # 带符号：>90° 表示核已经跑到雾的后面（连线反向）
    phi_signed = math.degrees(math.atan2(dSy, dSx)) if abs(dSx) > 1e-9 else 90.0
    # 观测口径：只能测「连线与碰撞轴的锐角」，不知道哪头是核
    phi_obs = phi_signed if phi_signed <= 90.0 else 180.0 - phi_signed
    return dict(dSx=dSx, dSy=dSy, dS=dS, phi=phi_obs, phi_signed=phi_signed,
                r=r * KMPS_MYR_TO_KPC, s=s * KMPS_MYR_TO_KPC, t_turn=t_turn)


def default_params():
    """一套全系统共用的参数（不是每个星系各调一套）。"""
    return dict(
        v0=2650.0,      # 初速度 km/s（Springel & Farrar 2007 量级）
        q=10.0,         # 主团/子团质量比
        tau=40.0,       # 雾的减速时间 Myr
        b=30.5,         # 回头前：碰撞参数留下的横向残余 kpc
                        # ★ 被 MACS(φ≤5°) 和 Bullet(φ≥11.8°) 夹在 30–31 的窄带里
        b_swing=230.0,  # 回头后：甩动幅度 kpc（saturate 模式用）
        tau_y=150.0,    # 甩动的时间尺度 Myr（saturate 模式用）
        swing_mode="power",  # power: 甩动随分离距离按幂律长，不饱和
        k_swing=3.2e-4,   # 幂律系数（按回程距离）
        alpha=2.28,     # 幂律指数（由 ZwCl + A520 两点定出）
        wd=1.0 / 220.0, # 振荡角频率（周期 ~1.4 Gyr）
        gamma=1.0 / 600.0,  # 阻尼（振幅衰减时间 600 Myr）
        brake="exp",
    )


def sample_curve(p, t_max=2000.0, n=4000):
    ts = [t_max * i / n for i in range(1, n + 1)]
    return [(t, trajectory(t, p)) for t in ts]


def report(p):
    out = []
    add = out.append
    add("=" * 74)
    add("四星系统一建模 — 《构造》核/雾相位轨迹")
    add("=" * 74)
    add("共用参数：v0=%.0f km/s  q=%.1f  tau=%.0f Myr" % (p["v0"], p["q"], p["tau"]))
    if p.get("swing_mode") == "power":
        add("横向两段式：回头前 b=%.1f kpc（线性长）→ 回头后甩动 "
            "dSy = b + %.2e·(回程距离)^%.2f（不饱和、连续）"
            % (p["b"], p["k_swing"], p["alpha"]))
    else:
        add("横向两段式：回头前 b=%.0f kpc（线性长）→ 回头后甩出 %.0f kpc，"
            "τy=%.0f Myr" % (p["b"], p.get("b_swing", 0.0), p["tau_y"]))
    add("相位刻度（非物理预测）：周期 = %.0f Myr   阻尼 = %.0f Myr"
        % (2 * math.pi / p["wd"], 1.0 / p["gamma"]))
    add("")

    # 核什么时候回头
    t_turn = math.pi / (2 * p["wd"])
    add("核第一次回头（速度过零）t = %.0f Myr" % t_turn)
    add("")

    add("-" * 74)
    add("相位采样：φ / ΔS 随时间的演化")
    add("-" * 74)
    add("   t(Myr) |   r(kpc) |   s(kpc) |   dSx   |   dSy   |   ΔS    "
        "|  φ观测  | 核位置")
    add("  ---------+----------+----------+---------+---------+---------"
        "+---------+--------")
    for t in (50, 100, 150, 200, 300, 400, 500, 700, 900, 1200, 1600):
        tr = trajectory(t, p)
        where = "雾前" if tr["dSx"] > 0 else "雾后"
        add("  %8.0f | %8.0f | %8.0f | %7.0f | %7.0f | %7.0f | %5.1f°  | %s"
            % (t, tr["r"], tr["s"], tr["dSx"], tr["dSy"], tr["dS"],
               tr["phi"], where))

    add("")
    add("-" * 74)
    add("从观测反推横向分量：dSy = ΔS · tan φ")
    add("-" * 74)
    add("  系统              ΔS中点   φ中点   反推 dSy   (kpc)")
    add("  -----------------+--------+--------+----------")
    for c in CLUSTERS:
        ds_mid = 0.5 * (c["dS_lo"] + c["dS_hi"])
        phi_mid = 0.5 * (c["phi_lo"] + c["phi_hi"])
        dsy = ds_mid * math.tan(math.radians(phi_mid))
        add("  %-17s %6.0f   %5.1f°   %8.0f"
            % (c["name"], ds_mid, phi_mid, dsy))
    add("")
    add("  ★ 两段式：MACS/Bullet 的 dSy 只有 20–30 kpc，")
    add("    A520/ZwCl 跳到 170–290 kpc —— 差一个数量级。")
    add("    这不像连续演化，更像跨过了一个**阈值**（核开始回头）。")

    add("")
    add("-" * 74)
    add("四个系统能否落在同一条曲线上")
    add("-" * 74)
    curve = sample_curve(p)
    for c in CLUSTERS:
        best = None
        for t, tr in curve:
            phi = tr["phi"]
            dS = tr["dS"]
            phi_ok = c["phi_lo"] <= phi <= c["phi_hi"]
            dS_ok = c["dS_lo"] <= dS <= c["dS_hi"]
            if phi_ok and dS_ok:
                best = (t, tr)
                break
        if best:
            t, tr = best
            add("  %-16s ✅ t=%4.0f Myr  φ=%5.1f°  ΔS=%4.0f kpc  (%s)"
                % (c["name"], t, tr["phi"], tr["dS"], c["phase"]))
        else:
            # 找最近的落点，报告差在哪
            near = min(curve, key=lambda kv: (
                max(0.0, kv[1]["phi"] - c["phi_hi"], c["phi_lo"] - kv[1]["phi"]) ** 2
                + (max(0.0, kv[1]["dS"] - c["dS_hi"], c["dS_lo"] - kv[1]["dS"]) / 10.0) ** 2))
            t, tr = near
            add("  %-16s ❌ 最近 t=%4.0f Myr  φ=%5.1f° (要 %.1f-%.1f)  "
                "ΔS=%4.0f kpc (要 %.0f-%.0f)  (%s)"
                % (c["name"], t, tr["phi"], c["phi_lo"], c["phi_hi"],
                   tr["dS"], c["dS_lo"], c["dS_hi"], c["phase"]))

    add("")
    add("-" * 74)
    add("压测：模型的自由度够不够")
    add("-" * 74)
    add("  φ 的两个来源：纵向被拉没（dSx 小） / 横向被甩出（dSy 大）")
    add("  共同点：两者都随时间**单调变大** —— 所以 φ 天然是时间的标签，")
    add("  四个系统可以按 φ 排序。真正的约束在 ΔS 上：")
    add("")
    add("  按 φ 排序：MACS(2-5) → Bullet(11.8) → A520(21.5-39) → ZwCl(32-41)")
    add("  对应 ΔS：  170-370   →   100-162   →    368-630     →   173-312")
    add("")
    add("  ΔS 不单调 —— 中间 A520 的分离距离最大。这是模型真正要回答的问题：")
    add("  如果 A520 已经在回程，为什么它反而是分得最开的？")
    add("")
    add("  两种可能的自洽读法：")
    add("    读法一：A520 不在回程，而是在「第二次穿出」的早期 —— 纵向重新拉开，")
    add("            但横向甩动已经积累（φ 大、ΔS 也大）。")
    add("    读法二：A520 的 ΔS 量法（368-630）混了中间介质/第三核，")
    add("            不是纯核-雾距离，不能和其他三个直接比。")
    add("  这两种读法给出不同的可检验预言：")
    add("    读法一 → A520 的雾应该在核**前面**（反向），或横向速度场有回摆特征。")
    add("    读法二 → 用统一的量法重测 A520，ΔS 应回落到 200-300 kpc。")
    return "\n".join(out)


def make_plot(p, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei",
                                   "PingFang SC", "Arial Unicode MS", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))

    # 左：φ–ΔS 平面，四个观测框 + 模型曲线
    ax = axes[0]
    curve = sample_curve(p, t_max=2000, n=600)
    phis = [tr["phi"] for _, tr in curve]
    dss = [tr["dS"] for _, tr in curve]
    ts = [t for t, _ in curve]
    sc = ax.scatter(dss, phis, c=ts, cmap="viridis", s=14, alpha=0.75, zorder=2)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("碰撞后时间 t (Myr)")

    for c in CLUSTERS:
        x0, x1 = c["dS_lo"], c["dS_hi"]
        y0, y1 = c["phi_lo"], c["phi_hi"]
        ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0,
                                   fill=True, alpha=0.22, color=c["color"],
                                   ec=c["color"], lw=1.5, zorder=3))
        ax.annotate(c["name"], xy=(0.5 * (x0 + x1), y1), xytext=(0, 5),
                    textcoords="offset points", ha="center", fontsize=9,
                    color=c["color"], fontweight="bold")

    ax.set_xlabel("核—雾分离距离 ΔS (kpc)")
    ax.set_ylabel("分离方向角 φ (°)")
    ax.set_title("四个并合系统 vs 《构造》相位轨迹（共用一套参数）")
    ax.grid(alpha=0.25)
    ax.set_xlim(0, 700)
    ax.set_ylim(0, 60)

    # 右：φ / 纵向 / 横向 随时间的演化
    ax = axes[1]
    ax.plot(ts, phis, lw=2.4, color="#c53030", label="φ (分离角)")
    ax.plot(ts, [tr["dSx"] for _, tr in curve], lw=1.8, ls="--",
            color="#2b6cb0", label="dSx 纵向 (kpc)")
    ax.plot(ts, [tr["dSy"] for _, tr in curve], lw=1.8, ls=":",
            color="#dd6b20", label="dSy 横向 (kpc)")
    ax.axhline(0, color="black", lw=0.8, alpha=0.4)

    for c in CLUSTERS:
        ymid = 0.5 * (c["phi_lo"] + c["phi_hi"])
        ax.axhspan(c["phi_lo"], c["phi_hi"], color=c["color"], alpha=0.10)
        ax.text(15, ymid, c["name"], fontsize=8, color=c["color"],
                va="center", fontweight="bold")

    ax.set_xlabel("碰撞后时间 t (Myr)")
    ax.set_ylabel("φ (°) / 分离分量 (kpc)")
    ax.set_title("φ 的两个来源：纵向被拉没、横向被甩出")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.25)

    fig.suptitle("《构造》四星系统一建模 —— 核/雾相位轨迹", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def main():
    ap = argparse.ArgumentParser(description="四星系统一建模")
    ap.add_argument("--plot", default=None, help="输出 PNG 路径")
    ap.add_argument("--out", default=None, help="报告写 UTF-8 文件")
    ap.add_argument("--v0", type=float, default=None)
    ap.add_argument("--q", type=float, default=None)
    ap.add_argument("--tau", type=float, default=None)
    ap.add_argument("--b", type=float, default=None)
    ap.add_argument("--tau-y", dest="tau_y", type=float, default=None)
    ap.add_argument("--period", type=float, default=None, help="振荡周期 Myr")
    ap.add_argument("--damp", type=float, default=None, help="阻尼时间 Myr")
    args = ap.parse_args()

    p = default_params()
    for k in ("v0", "q", "tau", "b", "tau_y"):
        if getattr(args, k) is not None:
            p[k] = getattr(args, k)
    if args.period:
        p["wd"] = 2 * math.pi / args.period
    if args.damp:
        p["gamma"] = 1.0 / args.damp

    text = report(p)
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
