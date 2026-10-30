#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
真实睡眠 EEG：频谱 vs 谐波（第二版，去伪迹）
================================================
数据：PhysioNet sleep-edf，SC4001E0-PSG.edf，EEG Fpz-Cz，100 Hz

第一版的问题：
  ① 2 秒窗 → 0.5 Hz 分辨率，分不清慢波和直流漂移
  ② 清醒段用了白天的活动段，眼动/肌肉伪迹造假双相干
  ③ 没高通滤波，1/f 背景压住一切

这一版的改进：
  ① 高通 0.5 Hz（Butterworth 4 阶，零相位）
  ② 4 秒窗（0.25 Hz 分辨率），50% 重叠
  ③ 清醒段自动挑 alpha（8–13 Hz）功率占比最高的 300 秒
  ④ 慢波段自动挑最长的连续 S3+S4
  ⑤ 伪迹剔除：分段峰峰值 > 150 μV 的段丢掉

判据不变：双相干 b(f) 高 → f 与 2f 相位锁定 → 真谐波。

用法：
    python eeg_bicoherence.py
    python eeg_bicoherence.py --plot ../figures/53_eeg_bicoherence.png
"""
import argparse
import os
import numpy as np
from scipy.signal import butter, filtfilt

_DATA = os.environ.get("EEG_DATA_DIR", ".")
PSG = os.path.join(_DATA, "SC4001E0-PSG.edf")
HYP = os.path.join(_DATA, "hyp.edf")
FS = 100.0
WIN = 4.0
OVERLAP = 0.5
HP = 0.5
ARTIFACT_UV = 150.0


def read_sig(start_s, dur_s, ch=0):
    import pyedflib
    f = pyedflib.EdfReader(PSG)
    sig = f.readSignal(ch, start=int(start_s * FS), n=int(dur_s * FS))
    f.close()
    return np.asarray(sig, dtype=float)


def highpass(x, fc=HP, fs=FS, order=4):
    b, a = butter(order, fc / (fs / 2), btype="high")
    return filtfilt(b, a, x)


def load_stages():
    import pyedflib
    f = pyedflib.EdfReader(HYP)
    ann = f.read_annotation()
    f.close()
    segs = []
    for a in ann:
        onset = a[0] / 1e7
        dur = float(a[1])
        label = a[2].decode()
        segs.append((onset, dur, label))
    return segs


def longest_run(segs, labels):
    best = (0.0, 0.0)
    i = 0
    while i < len(segs):
        j, total = i, 0.0
        while j < len(segs) and segs[j][2] in labels:
            total += segs[j][1]
            j += 1
        if total > best[0]:
            best = (total, segs[i][0])
        i = j + 1 if j == i else j
    return best


def alpha_ratio(x):
    n = int(4 * FS)
    if len(x) < n:
        return 0.0
    w = np.hanning(n)
    step = int(n * 0.5)
    tot, alpha = 0.0, 0.0
    for i in range(0, len(x) - n + 1, step):
        s = (x[i:i + n] - x[i:i + n].mean()) * w
        P = np.abs(np.fft.rfft(s)) ** 2
        f = np.fft.rfftfreq(n, 1 / FS)
        tot += P[(f > 1) & (f < 30)].sum()
        alpha += P[(f >= 8) & (f <= 13)].sum()
    return alpha / tot if tot > 0 else 0.0


def best_alpha_window(x, win_s=300.0, step_s=60.0):
    n = int(win_s * FS)
    step = int(step_s * FS)
    best = (-1.0, 0)
    for i in range(0, len(x) - n + 1, step):
        r = alpha_ratio(x[i:i + n])
        if r > best[0]:
            best = (r, i)
    return best[1] / FS, best[0]


def psd_welch(x, fs=FS, win=WIN):
    n = int(win * fs)
    step = int(n * (1 - OVERLAP))
    w = np.hanning(n)
    segs = []
    for i in range(0, len(x) - n + 1, step):
        s = x[i:i + n]
        if s.max() - s.min() > ARTIFACT_UV:
            continue
        s = s - s.mean()
        segs.append(np.fft.rfft(s * w))
    X = np.array(segs)
    P = np.mean(np.abs(X) ** 2, axis=0)
    freqs = np.fft.rfftfreq(n, 1 / fs)
    return freqs, P, X


def bicoherence_diag(X, freqs, fmax=20.0):
    idx = {round(f, 3): i for i, f in enumerate(freqs)}
    out_f, out_b = [], []
    for i, f in enumerate(freqs):
        if f < 1.0 or 2 * f > fmax:
            continue
        j = idx.get(round(2 * f, 3))
        if j is None:
            continue
        x1, x2 = X[:, i], X[:, j]
        B = np.mean(x1 * x1 * np.conj(x2))
        denom = np.sqrt(np.mean(np.abs(x1) ** 2) ** 2 * np.mean(np.abs(x2) ** 2))
        if denom <= 0:
            continue
        out_f.append(f)
        out_b.append(abs(B) / denom)
    return np.array(out_f), np.array(out_b)


def prepare():
    segs = load_stages()
    dur_slow, start_slow = longest_run(segs, {"Sleep stage 3", "Sleep stage 4"})
    x_slow = highpass(read_sig(start_slow, dur_slow))

    dur_w, start_w = longest_run(segs, {"Sleep stage W"})
    x_w_all = highpass(read_sig(start_w, dur_w))
    off, ratio = best_alpha_window(x_w_all, win_s=min(300.0, dur_w))
    x_w = x_w_all[int(off * FS):int(off * FS) + int(min(300.0, dur_w) * FS)]

    return {
        "慢波 S3+S4": (x_slow, dur_slow, start_slow),
        "清醒（alpha 最强）": (x_w, min(300.0, dur_w), start_w + off),
    }, ratio


def numbers():
    data, ar = prepare()
    print(f"清醒段 alpha 占比：{ar:.3f}")
    for name, (x, dur, start) in data.items():
        freqs, P, X = psd_welch(x)
        fb, b = bicoherence_diag(X, freqs)
        m = (freqs > 0.5) & (freqs < 20)
        peak = freqs[m][np.argmax(P[m])]
        print(f"{name}: {dur:.0f}s @ {start:.0f}s, 有效段 {len(X)}, "
              f"主峰 {peak:.2f} Hz, 双相干均值 {b.mean():.3f} (1–20 Hz)")


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    data, ar = prepare()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.4))

    for name, color in [("慢波 S3+S4", "#1f4e79"), ("清醒（alpha 最强）", "#c07000")]:
        x, dur, start = data[name]
        freqs, P, X = psd_welch(x)
        fb, b = bicoherence_diag(X, freqs)
        m = (freqs > 0.5) & (freqs < 30)
        ax1.semilogy(freqs[m], P[m], lw=1.6, color=color,
                     label=f"{name}（{dur:.0f}s, {len(X)} 段）")
        ax2.plot(fb, b, lw=1.8, color=color, label=name)

    ax1.set_xlabel("频率（Hz）")
    ax1.set_ylabel("功率（对数）")
    ax1.set_title("① 频谱（已高通 0.5 Hz）")
    ax1.legend(fontsize=9)
    ax1.grid(alpha=0.25, which="both")

    ax2.axhline(0.3, ls="--", color="#c00000", lw=1.0)
    ax2.text(1.1, 0.32, "随机信号显著线 ≈0.3", fontsize=8.5, color="#c00000")
    ax2.set_xlabel("频率 f（Hz），检测 f 与 2f 的相位锁定")
    ax2.set_ylabel("双相干 b(f)")
    ax2.set_title("② 有没有谐波结构")
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.25)
    ax2.set_ylim(0, 1)
    ax2.set_xlim(1, 20)

    fig.suptitle(f"真实睡眠 EEG（去伪迹版）：alpha 占比 {ar:.2f}", fontsize=13, y=1.0)
    fig.tight_layout()
    fig.savefig(path, dpi=145, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()
    if args.plot:
        os.makedirs(os.path.dirname(os.path.abspath(args.plot)), exist_ok=True)
        plot(args.plot)
        return
    numbers()


if __name__ == "__main__":
    main()
