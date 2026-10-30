#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
64 通道 EEG：alpha 是驻波还是行波？
========================================
数据：PhysioNet eegmmidb，S001R02.edf（闭眼静息，64 通道，160 Hz，61 s）

上一版的错误：单通道只能看"时间谐波"，看不到空间驻波。
驻波是空间模式 ψ(x,t)=A(x)cos(ωt)，A(x) 有节点和波腹。

判据：
    驻波  → 各通道相对某点的相位差集中在 0° 和 180° 两簇
    行波  → 相位差沿空间连续变化（0–360° 铺开）

做法：
    ① 去均值 + 共平均参考（CAR）
    ② 带通 8–13 Hz
    ③ Hilbert 取瞬时相位
    ④ 选 alpha 功率权重最高的时刻，算各通道相对于顶区(Cz)的相位差
    ⑤ 看分布：双峰还是铺开

★ 局限：只用了一个受试者、61 秒、一个频段；头模型和参考电极会影响相位。

用法：
    python eeg_standing_wave.py
    python eeg_standing_wave.py --plot ../figures/54_eeg_standing_wave.png
"""
import argparse
import os
import numpy as np
from scipy.signal import butter, filtfilt, hilbert

_DATA = os.environ.get("EEG_DATA_DIR", ".")
EDF = os.path.join(_DATA, "S001R02.edf")
BAND = (8.0, 13.0)


def load():
    import pyedflib
    f = pyedflib.EdfReader(EDF)
    fs = f.getSampleFrequency(0)
    labels = [f.getLabel(i).strip() for i in range(f.signals_in_file)]
    data = np.array([f.readSignal(i) for i in range(f.signals_in_file)])
    f.close()
    return data, labels, fs


def bandpass(x, band, fs, order=4):
    b, a = butter(order, [band[0] / (fs / 2), band[1] / (fs / 2)], btype="band")
    return filtfilt(b, a, x, axis=-1)


def analyse():
    data, labels, fs = load()
    data = data - data.mean(axis=1, keepdims=True)
    # 默认不做 CAR：共平均参考会强制相位归零，造出假的反相

    xb = bandpass(data, BAND, fs)
    amp = np.abs(hilbert(xb, axis=-1))
    ph = np.angle(hilbert(xb, axis=-1))

    # 取 alpha 功率最强的 10 秒
    w = int(10 * fs)
    env = amp.mean(axis=0)
    k = np.convolve(env, np.ones(w) / w, mode="same")
    i0 = int(np.argmax(k)) - w // 2
    i0 = max(0, min(i0, xb.shape[1] - w))
    sl = slice(i0, i0 + w)

    ref = labels.index("Cz..") if "Cz.." in labels else 0
    dphi = np.angle(np.exp(1j * (ph[:, sl] - ph[ref:ref + 1, sl])))
    wgt = amp[:, sl].mean(axis=1)
    mean_phase = np.angle(np.sum(wgt[:, None] * np.exp(1j * dphi), axis=1))
    return labels, mean_phase, wgt, fs


def numbers():
    labels, ph, wgt, fs = analyse()
    deg = np.rad2deg(ph)
    near0 = np.mean(np.abs(deg) < 45)
    near180 = np.mean(np.abs(np.abs(deg) - 180) < 45)
    mid = 1 - near0 - near180
    print(f"通道 {len(labels)}，alpha 带 {BAND[0]:.0f}–{BAND[1]:.0f} Hz")
    print(f"相对 Cz 的相位差分布：")
    print(f"  |Δφ|<45°（同相）      {near0*100:5.1f}%")
    print(f"  |Δφ|>135°（反相）     {near180*100:5.1f}%")
    print(f"  中间（90±45°）        {mid*100:5.1f}%")
    print("\n  驻波 → 前两项高、中间低；行波 → 三项均匀（各约 33%）")


def plot(path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams

    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    labels, ph, wgt, fs = analyse()
    deg = np.rad2deg(ph)

    fig = plt.figure(figsize=(13.5, 5.6))

    # ① 极坐标相位分布
    ax1 = fig.add_subplot(1, 2, 1, projection="polar")
    ax1.scatter(ph, wgt / wgt.max(), s=22, c="#1f4e79", alpha=0.7)
    ax1.set_title("① 各通道相对 Cz 的相位（半径＝alpha 幅度）", pad=18, fontsize=11)
    ax1.set_rticks([])

    # ② 直方图
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.hist(deg, bins=np.arange(-180, 181, 15), color="#1f4e79", alpha=0.75)
    ax2.axvline(0, ls="--", color="#888888", lw=1)
    ax2.axvline(180, ls="--", color="#888888", lw=1)
    ax2.axvline(-180, ls="--", color="#888888", lw=1)
    ax2.set_xlabel("相对 Cz 的相位差（度）")
    ax2.set_ylabel("通道数")
    ax2.set_title("② 相位差分布：双峰＝驻波，铺开＝行波", fontsize=11)
    ax2.grid(alpha=0.25)

    near0 = np.mean(np.abs(deg) < 45) * 100
    near180 = np.mean(np.abs(np.abs(deg) - 180) < 45) * 100
    mid = 100 - near0 - near180
    ax2.text(0.02, 0.95, f"同相 {near0:.0f}%\n反相 {near180:.0f}%\n中间 {mid:.0f}%",
             transform=ax2.transAxes, va="top", fontsize=9,
             bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    fig.suptitle("64 通道闭眼静息：alpha 是驻波还是行波？", fontsize=13, y=1.0)
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
