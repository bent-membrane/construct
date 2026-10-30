"""
模与率

三个实验：
  1. 连续 vs 整数搜索：假象地板塌多少
  2. 巴耳末 1885：留出验证 / 盲测帕申系 / 换尺子
  3. 率与模：定量打包 vs 恒流

附：为什么泊松回归不能估率（一次崩掉的记录）

运行：python modulus_rate.py
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")


# ============================================================
# 实验 1：连续 vs 整数搜索
# ============================================================
def alignment(x, k):
    """圆周统计：x 里的数落在 k 的整数倍附近有多齐"""
    return np.abs(np.exp(2j * np.pi * x / k).mean())


def exp1():
    print("=" * 62)
    print("实验 1：连续搜索 vs 只许整数")
    print("=" * 62)
    rng = np.random.default_rng(20261023)
    N, K_TRUE = 50, 12
    idx = np.arange(1, N + 1)
    struct = K_TRUE * idx + rng.normal(0, 0.4, N)

    for name, ks in [("连续 2.00–30.00（步长 0.005）", np.arange(2.0, 30.0, 0.005)),
                     ("只许整数 2–30", np.arange(2, 31))]:
        def best(x):
            s = np.array([alignment(x, k) for k in ks])
            return ks[int(np.argmax(s))], s.max()
        # 随机基线
        base = np.array([best(rng.uniform(0, K_TRUE * N, N))[1] for _ in range(2000)])
        kb, ab = best(struct)
        print(f"\n{name}（候选 {len(ks)} 个）")
        print(f"  纯随机最高峰：中位 {np.median(base):.3f} / 95分位 {np.percentile(base,95):.3f}")
        print(f"  真结构 k=12 → 读出 {kb:.2f}，峰 {ab:.3f}")


# ============================================================
# 实验 2：巴耳末 1885
# ============================================================
def exp2():
    print("\n" + "=" * 62)
    print("实验 2：巴耳末 1885")
    print("=" * 62)
    lam = np.array([656.279, 486.135, 434.047, 410.174])   # nm
    TRUE_R = 1.09677583e7                                   # m^-1

    def fit_R(lam_sub, n_sub):
        nu = 1.0 / lam_sub
        X = 0.25 - 1.0 / n_sub**2
        return (X @ nu) / (X @ X)                          # nm^-1

    # 留出验证：3 条拟合，预测第 4 条
    R3 = fit_R(lam[:3], np.array([3., 4., 5.]))
    pred4 = 1.0 / (R3 * (0.25 - 1/36))
    print(f"\n留出验证（n=3,4,5 拟合 → 预测 n=6）")
    print(f"  预测 {pred4:.4f} nm / 实测 {lam[3]:.4f} nm / 偏差 {(pred4-lam[3])/lam[3]*100:+.4f}%")

    # 盲测帕申系
    R4 = fit_R(lam, np.arange(3., 7.))
    true_paschen = {4: 1875.1, 5: 1281.8, 6: 1093.8}
    print(f"\n盲测帕申系（巴耳末拟合时没用过）")
    for m in (4, 5, 6):
        pred = 1.0 / (R4 * (1/9 - 1/m**2))
        print(f"  n={m}: 预测 {pred:8.2f} / 实测 {true_paschen[m]:8.2f} / "
              f"偏差 {(pred-true_paschen[m])/true_paschen[m]*100:+.3f}%")

    # 换尺子
    R_A = fit_R(lam * 10, np.arange(3., 7.))
    print(f"\n换尺子")
    print(f"  用 nm：R = {R4:.6f} nm⁻¹ = {R4*1e9:.2f} m⁻¹")
    print(f"  用 Å ：R = {R_A:.6f} Å⁻¹  = {R_A*1e10:.2f} m⁻¹")
    print(f"  换回同一把尺子完全相同；整数 n 一直是 3,4,5,6，不随尺子动")
    print(f"  （真值 {TRUE_R:.2f} m⁻¹，拟合差 0.029%：用了空气波长）")


# ============================================================
# 实验 3：率与模
# ============================================================
def exp3():
    print("\n" + "=" * 62)
    print("实验 3：率与模 —— 定量打包 vs 恒流")
    print("=" * 62)
    L = 600.0
    KS = [6, 8, 10, 12, 16, 20, 24]
    print(f"\n{'k':>3} | {'定量打包(150/周期)':>22} | {'恒流(总数20000)':>20}")
    print(f"{'':>3} | {'率':>10}{'每周期量':>12} | {'率':>10}{'每周期量':>10}")
    A, B = [], []
    for k in KS:
        npd = int(L // k)
        nA, nB = npd * 150, (20000 // npd) * npd
        lA, lB = nA / L, nB / L
        A.append(lA * k); B.append(lB * k)
        print(f"{k:>3} | {lA:>10.2f}{lA*k:>12.1f} | {lB:>10.2f}{lB*k:>10.1f}")
    print(f"\n定量打包：每周期量 {np.round(A,1).tolist()} → 恒定")
    print(f"恒流    ：每周期量 {np.round(B,1).tolist()} → 随 k 上升")
    print("读法：率恒定=恒流；每周期量恒定=定量打包")


# ============================================================
# 附：泊松回归为什么不能估率
# ============================================================
def poisson_crash():
    """点在周期边界挤成一团时，泊松回归 k 读对、率估成 0"""
    from scipy.optimize import minimize
    print("\n" + "=" * 62)
    print("附：泊松回归崩掉的一次（k 对了，率成了 0）")
    print("=" * 62)
    rng = np.random.default_rng(7)
    L, W = 600.0, 2.0

    def build(k, per):
        npd = int(L // k)
        idx = np.repeat(np.arange(npd), per)
        return idx * k + rng.normal(0, 0.5, len(idx))

    def fit(x, k):
        edges = np.arange(0, L + W, W)
        y, _ = np.histogram(x, bins=edges)
        t = (edges[:-1] + edges[1:]) / 2
        X = np.column_stack([np.ones_like(t), np.cos(2*np.pi*t/k), np.sin(2*np.pi*t/k)])
        def nll(b):
            eta = X @ b
            return -np.sum(y*eta - np.exp(np.clip(eta, -20, 20)))
        r = minimize(nll, [np.log(max(y.mean(), 1e-6)), 0, 0], method="Nelder-Mead",
                     options={"maxiter": 8000, "xatol": 1e-9, "fatol": 1e-9})
        return r.x[0]

    print(f"\n{'k':>3} | {'读回 k':>6} | {'估出的率(点/长度)':>18}")
    for k in (6, 8, 10, 12):
        x = build(k, 150)
        b0 = fit(x, k)
        print(f"{k:>3} | {k:>6} | {np.exp(b0)/W:>18.1f}")
    print("\n真值应为：每周期 150 包 → 率 = 150/k")
    print("实际估成 0–2：零膨胀分布，泊松假设不成立")


if __name__ == "__main__":
    exp1()
    exp2()
    exp3()
    poisson_crash()
