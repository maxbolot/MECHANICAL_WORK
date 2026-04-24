import matplotlib.pyplot as plt
import numpy as np


def plot_percentile_bin_ratio_comparison(
    stats_ctrl,
    stats_warm,
    region_title: str,
    label_ctrl: str = "Control",
    label_warm: str = "+4K",
):
    x = np.arange(1, len(stats_ctrl.ratio_bin) + 1)
    y_ctrl = np.asarray(stats_ctrl.ratio_bin, dtype=float)
    y_warm = np.asarray(stats_warm.ratio_bin, dtype=float)

    color_ctrl = (0.20, 0.45, 0.70)
    color_warm = (0.80, 0.20, 0.20)
    edge_ctrl = (0.05, 0.15, 0.25)
    edge_warm = (0.35, 0.05, 0.05)

    fig, ax = plt.subplots(figsize=(12, 5), dpi=110)

    for ib, xpos in enumerate(x):
        if np.isfinite(y_warm[ib]):
            ax.bar(xpos, y_warm[ib], width=0.75, color=color_warm, edgecolor=edge_warm, alpha=1.0)
        if np.isfinite(y_ctrl[ib]):
            ax.bar(xpos, y_ctrl[ib], width=0.75, color=color_ctrl, edgecolor=edge_ctrl, alpha=0.75)

    (h_ctrl,) = ax.plot(x, y_ctrl, ".-", color=tuple(v * 0.7 for v in color_ctrl), linewidth=1.2, markersize=12)
    (h_warm,) = ax.plot(x, y_warm, ".-", color=tuple(v * 0.7 for v in color_warm), linewidth=1.2, markersize=12)

    ax.set_xticks(x)
    ax.set_xticklabels(stats_ctrl.labels, rotation=35, ha="right")
    ax.set_xlim(0.5, len(x) + 0.5)
    ax.set_xlabel("Percentile-rank bin")
    ax.set_ylabel("Lift/Work ratio")
    ax.set_title(f"{region_title}: Lift/Work Ratio by Percentile-Rank Bin ({label_ctrl} vs {label_warm})")
    ax.grid(True)
    ax.legend([h_ctrl, h_warm], [label_ctrl, label_warm], loc="upper left")

    fig.tight_layout()
    return fig, ax
