#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
子弹星团沙盘 —— 核/雾位移比与分离方向
=======================================

把《构造》里"核=断点、雾=粉尘"的碰撞图像落到一个可算的玩具模型上。

核心判据（不依赖碰撞时间）：
    ΔS / S_core = 1 - <v_fog> / v_core

真实子弹星团 (1E 0657-56) 的靶值（个案，只标量级）：
    ΔS/S_core   = 0.25 - 0.60          （随距离量法 100/162 kpc 与碰撞后时间变）
    <v_fog>/v_core = 0.40 - 0.75       （由上式的平均速度比，不是瞬时速度）

物理设定：
    - 核（断点）无阻力，匀速 v_core，沿碰撞轴运动
    - 雾（气体）从 v_core 减速到质心速度 v_inf = v_core / q
      （q = 主团/子团质量比，质心速度在主团系里约为 v_core/q）
    - 阻力模型决定 v(t) 的形状：
        exp   指数衰减   dv/dt = -(v - v_inf)/tau
        quad  平方阻力   dv/dt = -(v - v_inf)^2 / ((v0 - v_inf) tau)
        const 匀减速     tau = 从 v0 降到 v_inf 所需时间
    - 位移比只依赖 tau/t，不依赖 t 的绝对值

脱落比例（--f）：
    - 尘分两部分：f 比例被震落、脱离后减速到 v_inf；1-f 仍附着，跟着核走
    - ΔS/S_core = f * dr0，dr0 = 全部脱落（f=1）时的位移比
    - 观测 0.25-0.60 反推：若脱落部分完全减速（dr0_max = 1-1/q），
      f >= 0.28-0.67；减速不完全时 dr0 更小、f 更高，最高到 1
    - 标准模型没有这个自由度（气体本来就不附着，f=1）

碰撞参数与分离方向（--b / --eta / --inc）：
    - b   两团轨迹的碰撞参数，kpc（横向错开的距离）
    - eta 雾的横向耦合系数：雾在 t 时间内被侧推 eta*b（0 = 纯纵向分离）
    - inc 轨道相对天空平面的倾角；纵向分离投影缩短 cos(inc)，横向不变
    - 分离方向角  phi = atan( eta*b / (ΔS_x * cos inc) )，相对碰撞轴
    - 观测量是 phi；子弹星团的核-雾连线接近沿碰撞轴，所以 phi_obs 很小，
      这反过来给 eta 一个上限。引力偏折忽略（核轨迹近似直线）。

用法：
    python bullet_cluster_sandbox.py
    python bullet_cluster_sandbox.py --q 6 --model quad
    python bullet_cluster_sandbox.py --tau 40
    python bullet_cluster_sandbox.py --b 100 --eta 0.3 --inc 10
    python bullet_cluster_sandbox.py --f 0.75
    python bullet_cluster_sandbox.py --out report.txt   # 报告写 UTF-8 文件
    python bullet_cluster_sandbox.py --projection proj.png   # 3D 投影：φ_obs 随视线倾角 i
"""

import argparse
import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# 1 km/s × 1 Myr = 3.15576e13 km = 1.0227e-3 kpc
KMPS_MYR_TO_KPC = 1.0227e-3

MODELS = ("exp", "quad", "const")

# 个案参考带（子弹星团；只标量级，不钉死）
TARGET_LO = 0.25
TARGET_HI = 0.60
PHI_OBS = 11.0   # 个案 3D 角（Clowe 2006 坐标算出 11.8-14.1°，去投影后约 11°）


def fog_v(v0, vinf, t, tau, model):
    """雾在 t 时刻的速度 (km/s)。"""
    if model == "exp":
        return vinf + (v0 - vinf) * math.exp(-t / tau)
    if model == "quad":
        return vinf + (v0 - vinf) / (1.0 + t / tau)
    if model == "const":
        if tau <= 0 or t >= tau:
            return vinf
        return v0 - (v0 - vinf) * t / tau
    raise ValueError("unknown model: %s" % model)


def fog_s(v0, vinf, t, tau, model, n=4000):
    """雾在 t 时间内的位移 (km/s · Myr)，中点法数值积分。"""
    dt = t / n
    total = 0.0
    for i in range(n):
        total += fog_v(v0, vinf, (i + 0.5) * dt, tau, model)
    return total * dt


def evaluate(vcore, t, q, model, tau, b=0.0, eta=0.0, inc=0.0, f=1.0):
    """返回位移比、瞬时速度比、分离分量和方向角。f = 震落比例。"""
    vinf = vcore / q
    vfog_det = fog_v(vcore, vinf, t, tau, model)   # 震落部分的瞬时速度
    sfog_det = fog_s(vcore, vinf, t, tau, model)   # 震落部分的位移
    score = vcore * t

    # 混合：f 震落（被拖慢），1-f 仍附着（跟着核走）
    vfog = f * vfog_det + (1.0 - f) * vcore
    sfog = f * sfog_det + (1.0 - f) * score

    vr = vfog / vcore
    dr = 1.0 - sfog / score
    dSx = (score - sfog) * KMPS_MYR_TO_KPC          # 纵向分离, kpc
    dSy = f * eta * b                               # 横向分离, kpc
    dSx_proj = dSx * math.cos(math.radians(inc))    # 天空平面上的纵向分量
    if dSx_proj > 0:
        phi = math.degrees(math.atan2(dSy, dSx_proj))
    else:
        phi = 90.0 if dSy else 0.0

    return {
        "vr": vr,
        "dr": dr,
        "vfog": vfog,
        "score_kpc": score * KMPS_MYR_TO_KPC,
        "sfog_kpc": sfog * KMPS_MYR_TO_KPC,
        "dSx": dSx,
        "dSy": dSy,
        "phi": phi,
    }


def solve_tau(target_ratio, vcore, t, q, model):
    """二分求使 ΔS/S_core = target_ratio 的 tau（比值随 tau 单调递减）。"""
    lo, hi = 1e-3, 1e5
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if evaluate(vcore, t, q, model, mid)["dr"] > target_ratio:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def solve_tau_for_speed(target_vr, vcore, t, q, model):
    """二分求使瞬时 v_fog/v_core = target_vr 的 tau（随 tau 单调递增）。"""
    lo, hi = 1e-3, 1e5
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if evaluate(vcore, t, q, model, mid)["vr"] < target_vr:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def build_report(args):
    lines = []
    add = lines.append

    vcore, t, q = args.vcore, args.t, args.q
    b, eta, inc = args.b, args.eta, args.inc
    vinf = vcore / q
    add("=" * 72)
    add("Bullet Cluster sandbox -- core/fog displacement + separation angle")
    add("=" * 72)
    add("v_core = %.0f km/s   t = %.0f Myr   q = %.1f   v_inf = %.0f km/s"
        % (vcore, t, q, vinf))
    add("S_core = %.0f kpc   b = %.0f kpc   eta = %.2f   inc = %.1f deg   f = %.2f"
        % (vcore * t * KMPS_MYR_TO_KPC, b, eta, inc, args.f))
    add("")
    add("target band:  dS/S_core in [%.2f, %.2f]   mean v_fog/v_core in [%.2f, %.2f]"
        % (args.target_lo, args.target_hi, 1.0 - args.target_hi, 1.0 - args.target_lo))
    add("")

    models = MODELS if args.model == "all" else (args.model,)

    if args.tau is not None:
        add("-" * 72)
        add("single point: tau = %.1f Myr" % args.tau)
        add("-" * 72)
        for m in models:
            r = evaluate(vcore, t, q, m, args.tau, b, eta, inc, args.f)
            add("%-6s  v/v = %.3f (%6.0f km/s)  dS/S = %.3f  "
                "dSx = %5.0f kpc  phi = %5.1f deg"
                % (m, r["vr"], r["vfog"], r["dr"], r["dSx"], r["phi"]))
        return "\n".join(lines)

    # 扫描
    add("-" * 72)
    add("scan: tau (Myr) -> [instantaneous v_fog/v_core, dS/S_core]")
    add("-" * 72)
    add("  tau |" + "".join(" %-14s" % m for m in models))
    add("  ----+" + "-" * (15 * len(models)))
    for tau in (20, 30, 40, 50, 75, 100, 150, 200, 300):
        row = " %4d |" % tau
        for m in models:
            r = evaluate(vcore, t, q, m, tau, f=args.f)
            row += " %5.2f / %-6.2f" % (r["vr"], r["dr"])
        add(row)

    # 反演：靶值对应的 tau
    add("")
    add("-" * 72)
    add("required tau (Myr) to hit the target band")
    add("-" * 72)
    vr_hi = 1.0 - args.target_lo
    vr_lo = 1.0 - args.target_hi
    add("  model | tau for dS/S=%.2f | tau for dS/S=%.2f | tau for v/v=%.2f | tau for v/v=%.2f"
        % (args.target_hi, args.target_lo, vr_hi, vr_lo))
    add("  ------+---------------------+--------------------+-------------------+------------------")
    for m in models:
        tau_hi = solve_tau(args.target_hi, vcore, t, q, m)
        tau_lo = solve_tau(args.target_lo, vcore, t, q, m)
        tau_v_hi = solve_tau_for_speed(vr_hi, vcore, t, q, m)
        tau_v_lo = solve_tau_for_speed(vr_lo, vcore, t, q, m)
        add("  %-6s| %19.1f | %18.1f | %17.1f | %16.1f"
            % (m, tau_hi, tau_lo, tau_v_hi, tau_v_lo))

    # 分离方向角
    add("")
    add("-" * 72)
    add("separation angle on the sky: phi = atan( eta*b / (dSx*cos(inc)) )")
    add("-" * 72)
    add("  tau |" + "".join(" %-12s" % m for m in models))
    add("  ----+" + "-" * (13 * len(models)))
    for tau in (20, 30, 50, 75, 100, 150):
        row = " %4d |" % tau
        for m in models:
            r = evaluate(vcore, t, q, m, tau, b, eta, inc, args.f)
            row += " %7.1f deg" % r["phi"]
        add(row)

    # 方向约束：phi_obs < phimax 时 eta 的上限
    add("")
    add("-" * 72)
    add("direction constraint: phi_obs < %.0f deg  ->  eta_max = tan(phi)*dSx/b"
        % args.phimax)
    add("-" * 72)
    for m in models:
        tau_mid = solve_tau(0.5 * (args.target_lo + args.target_hi), vcore, t, q, m)
        r = evaluate(vcore, t, q, m, tau_mid, b, eta, inc, args.f)
        if b > 0:
            eta_max = math.tan(math.radians(args.phimax)) * r["dSx"] / b
            add("  %-6s tau=%6.1f Myr  dSx=%5.0f kpc  ->  eta_max = %.2f"
                % (m, tau_mid, r["dSx"], eta_max))
        else:
            add("  %-6s tau=%6.1f Myr  dSx=%5.0f kpc  ->  eta_max = inf (b=0)"
                % (m, tau_mid, r["dSx"]))

    # 脱落比例
    add("")
    add("-" * 72)
    add("detachment fraction f:  dS/S_core = f * dr0")
    add("-" * 72)
    add("  dr0 = displacement ratio if ALL dust detaches (f=1)")
    add("  fully-braked limit: dr0_max = 1 - 1/q = %.2f" % (1.0 - 1.0 / q))
    add("  observed dS/S = %.2f-%.2f  ->  f >= %.2f - %.2f  (lower bound)"
        % (args.target_lo, args.target_hi,
           args.target_lo / (1.0 - 1.0 / q),
           args.target_hi / (1.0 - 1.0 / q)))
    add("  incomplete braking gives smaller dr0 and larger f, up to f=1")
    add("")
    add("  tau |" + "".join(" %-15s" % m for m in models))
    add("  ----+" + "-" * (16 * len(models)))
    for tau in (20, 30, 50, 75, 100, 150):
        row = " %4d |" % tau
        for m in models:
            dr0 = evaluate(vcore, t, q, m, tau)["dr"]
            fneed = 0.5 * (args.target_lo + args.target_hi) / dr0 if dr0 > 0 else float("inf")
            row += " dr0=%.2f f=%5.2f" % (dr0, fneed)
        add(row)
    add("  (f>1 = even full detachment cannot reach the band midpoint at that tau)")
    add("")
    add("  note: phi is independent of f -- f scales dSx and dSy together and")
    add("  cancels in atan(dSy/dSx). Displacement ratio measures f; the")
    add("  separation angle measures eta*b. The two observables are orthogonal.")

    add("")
    add("interpretation:")
    add("  - v_fog/v_core above is the INSTANTANEOUS speed at t; dS/S_core is the")
    add("    time-averaged displacement ratio. The observable is the latter.")
    add("  - reproducing dS/S = %.2f-%.2f requires:" % (args.target_lo, args.target_hi))
    for m in models:
        t1 = solve_tau(args.target_lo, vcore, t, q, m)
        t2 = solve_tau(args.target_hi, vcore, t, q, m)
        add("      %-5s tau ~ %.0f-%.0f Myr" % (m, min(t1, t2), max(t1, t2)))
    add("  - only quadratic (ram-pressure) braking leaves a residual speed near")
    add("    the observed mean, consistent with still-moving X-ray gas.")
    add("  - phi is a SECOND observable: the Bullet Cluster's DM-gas separation")
    add("    lies close to the collision axis, so phi_obs is small and eta*b must")
    add("    stay well below dSx (see eta_max above).")
    return "\n".join(lines)


def make_plot(args, path):
    """四联图：轨迹、位移比-τ、方向角-τ、位移比-f。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = [
        "Microsoft YaHei", "SimHei", "PingFang SC", "Heiti SC",
        "Arial Unicode MS", "DejaVu Sans",
    ]
    rcParams["axes.unicode_minus"] = False

    vcore, t_end, q = args.vcore, args.t, args.q
    b, eta, inc, f = args.b, args.eta, args.inc, args.f
    tau_ref = args.tau if args.tau is not None else 40.0
    model_ref = args.model if args.model != "all" else "exp"

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # 1) 轨迹
    ax = axes[0][0]
    n = 400
    ts = [t_end * i / n for i in range(n + 1)]
    x_core = [vcore * tt * KMPS_MYR_TO_KPC for tt in ts]
    y_core = [b for _ in ts]
    x_fog, y_fog = [], []
    for tt in ts:
        sd = fog_s(vcore, vcore / q, tt, tau_ref, model_ref)
        sf = f * sd + (1.0 - f) * vcore * tt
        x_fog.append(sf * KMPS_MYR_TO_KPC)
        y_fog.append(b + f * eta * b * (tt / t_end))
    ax.plot(x_core, y_core, "-", lw=2.2, color="#2b6cb0", label="核（断点，无阻力）")
    ax.plot(x_fog, y_fog, "-", lw=2.2, color="#dd6b20", label="雾（尘，部分震落）")
    ax.axhline(b, color="#2b6cb0", ls=":", lw=1, alpha=0.5)
    ax.plot([0], [b], "o", color="black", ms=6)
    ax.annotate("碰撞点", (0, b), textcoords="offset points", xytext=(8, 8))
    ax.text(0.98, 0.95, "注意：x/y 轴比例不同\n图中角度被放大，模型 φ = %.1f°"
            % evaluate(vcore, t_end, q, model_ref, tau_ref, b, eta, inc, f)["phi"],
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            color="#718096")
    ax.annotate("", xy=(x_core[-1], y_core[-1]), xytext=(x_fog[-1], y_fog[-1]),
                arrowprops=dict(arrowstyle="<->", color="gray", lw=1.6))
    r = evaluate(vcore, t_end, q, model_ref, tau_ref, b, eta, inc, f)
    ax.text(0.02, 0.97,
            "ΔSx = %.0f kpc\nΔSy = %.0f kpc\nφ = %.1f°（模型）"
            % (r["dSx"], r["dSy"], r["phi"]),
            transform=ax.transAxes, va="top", fontsize=11,
            bbox=dict(boxstyle="round", fc="white", ec="#cccccc", alpha=0.9))
    ax.set_xlabel("沿碰撞轴 x (kpc)")
    ax.set_ylabel("横向 y (kpc)")
    ax.set_title("轨迹（t = %.0f Myr, %s, τ = %.0f Myr）" % (t_end, model_ref, tau_ref))
    ax.legend(loc="center right")
    ax.grid(alpha=0.25)

    # 2) 位移比 vs τ
    ax = axes[0][1]
    taus = [i * 2.0 for i in range(2, 126)]
    for m in MODELS:
        vals = [evaluate(vcore, t_end, q, m, tau, f=f)["dr"] for tau in taus]
        ax.plot(taus, vals, lw=2, label=m)
    ax.axhspan(args.target_lo, args.target_hi, color="#38a169", alpha=0.15, label="观测参考带")
    ax.set_xlabel("减速时间 τ (Myr)")
    ax.set_ylabel("ΔS / S_core")
    ax.set_title("位移比 vs 减速时间（f = %.2f）" % f)
    ax.legend()
    ax.grid(alpha=0.25)

    # 3) 方向角 vs τ
    ax = axes[1][0]
    for m in MODELS:
        vals = [evaluate(vcore, t_end, q, m, tau, b, eta, inc, f)["phi"] for tau in taus]
        ax.plot(taus, vals, lw=2, label=m)
    ax.axhspan(0, args.phimax, color="#38a169", alpha=0.12,
               label="φ < %.0f° 观测上限" % args.phimax)
    ax.axhline(PHI_OBS, color="#c53030", ls="-.", lw=1.8,
               label="个案 3D φ = %.1f°（去投影）" % PHI_OBS)
    ax.set_xlabel("减速时间 τ (Myr)")
    ax.set_ylabel("分离方向角 φ (°)")
    ax.set_title("方向角 vs 减速时间（b = %.0f kpc, η = %.2f）" % (b, eta))
    ax.legend()
    ax.grid(alpha=0.25)

    # 4) 位移比 vs f
    ax = axes[1][1]
    fs = [i / 50.0 for i in range(51)]
    for tau0 in (20, 50, 100):
        vals = [evaluate(vcore, t_end, q, "exp", tau0, f=ff)["dr"] for ff in fs]
        ax.plot(fs, vals, lw=2, label="τ = %.0f Myr" % tau0)
    ax.axhspan(args.target_lo, args.target_hi, color="#38a169", alpha=0.15, label="观测参考带")
    ax.set_xlabel("震落比例 f")
    ax.set_ylabel("ΔS / S_core")
    ax.set_title("位移比 vs 震落比例（exp 模型）")
    ax.legend()
    ax.grid(alpha=0.25)

    fig.suptitle("子弹星团沙盘 — 《构造》核/雾判据", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("plot written: %s" % path)


def model_phi_true(args):
    """当前参数下的真实分离角（不含投影），度。"""
    model_ref = args.model if args.model != "all" else "exp"
    tau_ref = args.tau if args.tau is not None else 40.0
    r = evaluate(args.vcore, args.t, args.q, model_ref, tau_ref,
                 args.b, args.eta, 0.0, args.f)
    if r["dSx"] <= 0:
        return 0.0
    return math.degrees(math.atan2(r["dSy"], r["dSx"]))


def phi_obs_from_i(phi_true_deg, i_deg):
    """真实角按视线倾角 i 投影到天空平面。i=0 侧视（不压缩），i=90 俯视（纵向压没）。"""
    c = math.cos(math.radians(i_deg))
    if c <= 1e-9:
        return 90.0
    return math.degrees(math.atan(math.tan(math.radians(phi_true_deg)) / c))


def print_projection_table(args):
    phi_model = model_phi_true(args)
    print("model phi (no projection)       = %.2f deg" % phi_model)
    print("case 3D phi (deprojected)       = %.2f deg" % PHI_OBS)
    print("%8s %10s %10s %10s" % ("i(deg)", "phi_model", "phi_case", "dSx kept"))
    for i in (0, 10, 20, 30, 35, 45, 60, 65, 75, 90):
        kept = "0%" if i == 90 else "%.0f%%" % (100 * math.cos(math.radians(i)))
        print("%8.0f %10.2f %10.2f %10s"
              % (i, phi_obs_from_i(phi_model, i),
                 phi_obs_from_i(PHI_OBS, i), kept))


def make_projection_plot(args, path):
    """3D 投影：视线倾角 i 如何扭曲观测到的分离角。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = [
        "Microsoft YaHei", "SimHei", "PingFang SC", "Heiti SC",
        "Arial Unicode MS", "DejaVu Sans",
    ]
    rcParams["axes.unicode_minus"] = False

    phi_true = model_phi_true(args)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))
    ii = [k * 0.5 for k in range(181)]

    # 左：曲线族
    ax = axes[0]
    for pt in (3.0, 6.3, 10.0, 15.0):
        ax.plot(ii, [phi_obs_from_i(pt, x) for x in ii], lw=1.8,
                label="真实 φ = %.1f°" % pt)
    ax.plot(ii, [phi_obs_from_i(PHI_OBS, x) for x in ii], lw=2.8, ls="-.",
            color="#c53030", label="个案 3D φ = %.1f°（去投影）" % PHI_OBS)
    ax.plot(ii, [phi_obs_from_i(phi_true, x) for x in ii], lw=3.2, ls="--",
            color="black", label="当前模型 φ = %.1f°" % phi_true)
    ax.axvspan(10, 20, color="#3182ce", alpha=0.16, label="文献倾角 10–20°（68%）")
    ax.axvspan(0, 35, color="#3182ce", alpha=0.06, label="95% 上限 35°")
    ax.axhline(15, color="#dd6b20", ls=":", lw=1.6, label="15° 观测上限")
    ax.set_xlabel("视线倾角 i (°)   —   0 = 侧视，90 = 俯视")
    ax.set_ylabel("观测到的夹角 φ_obs (°)")
    ax.set_title("投影对可见角的扭曲")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.25)

    # 右：散点
    ax = axes[1]
    pts_i = [k * 2 for k in range(46)]
    colors = {3.0: "#718096", 6.3: "#3182ce", 10.0: "#38a169", 15.0: "#dd6b20"}
    for pt in (3.0, 6.3, 10.0, 15.0):
        ax.scatter(pts_i, [phi_obs_from_i(pt, x) for x in pts_i],
                   s=22, alpha=0.75, color=colors[pt], label="真实 φ = %.1f°" % pt)
    ax.scatter(pts_i, [phi_obs_from_i(PHI_OBS, x) for x in pts_i],
               s=26, alpha=0.85, color="#c53030",
               label="个案 3D φ = %.1f°（去投影）" % PHI_OBS)
    ax.axvspan(10, 20, color="#3182ce", alpha=0.16)
    ax.axhline(15, color="#dd6b20", ls=":", lw=1.6)
    if PHI_OBS < 15:
        need_i = math.degrees(math.acos(
            math.tan(math.radians(PHI_OBS)) / math.tan(math.radians(15))))
        ax.axvline(need_i, color="black", ls="--", lw=1.2)
        ax.annotate("个案 3D φ=%.1f° → 15° 需 i=%.0f°" % (PHI_OBS, need_i),
                    xy=(need_i, 15),
                    xytext=(max(need_i - 34, 2), 8),
                    fontsize=9,
                    arrowprops=dict(arrowstyle="->", color="black", lw=1))
    ax.set_xlabel("视线倾角 i (°)")
    ax.set_ylabel("φ_obs (°)")
    ax.set_title("散点：i 每 2° 采样一次")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.25)

    fig.suptitle("3D 投影：碰撞轴倾斜如何扭曲核—雾分离角", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("projection plot written: %s" % path)


def main():
    p = argparse.ArgumentParser(description="Bullet Cluster core/fog sandbox")
    p.add_argument("--vcore", type=float, default=2650.0,
                   help="core speed km/s (default 2650, Springel & Farrar 2007)")
    p.add_argument("--t", type=float, default=150.0,
                   help="time since core passage, Myr (default 150, Markevitch)")
    p.add_argument("--q", type=float, default=10.0,
                   help="mass ratio main/sub (default 10)")
    p.add_argument("--model", default="all", choices=("exp", "quad", "const", "all"))
    p.add_argument("--tau", type=float, default=None,
                   help="fog braking time in Myr; omit to scan")
    p.add_argument("--b", type=float, default=100.0,
                   help="impact parameter in kpc (default 100)")
    p.add_argument("--eta", type=float, default=0.3,
                   help="fog transverse coupling: dy = eta*b (default 0.3)")
    p.add_argument("--inc", type=float, default=0.0,
                   help="orbit inclination to sky plane, deg (default 0)")
    p.add_argument("--f", type=float, default=1.0,
                   help="fraction of dust knocked off the membrane (default 1.0)")
    p.add_argument("--phimax", type=float, default=15.0,
                   help="observed upper limit on separation angle, deg (default 15)")
    p.add_argument("--target-lo", type=float, default=TARGET_LO,
                   help="lower edge of dS/S_core reference band (default 0.25)")
    p.add_argument("--target-hi", type=float, default=TARGET_HI,
                   help="upper edge of dS/S_core reference band (default 0.60)")
    p.add_argument("--out", default=None, help="write report to UTF-8 file")
    p.add_argument("--plot", default=None, help="write a 4-panel PNG to this path")
    p.add_argument("--projection", default=None,
                   help="write a 3D-projection PNG (phi_obs vs viewing angle i)")
    args = p.parse_args()

    report = build_report(args)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(report + "\n")
        print("report written: %s" % args.out)
    elif not args.plot and not args.projection:
        print(report)

    if args.plot:
        make_plot(args, args.plot)

    if args.projection:
        print_projection_table(args)
        make_projection_plot(args, args.projection)


if __name__ == "__main__":
    main()
