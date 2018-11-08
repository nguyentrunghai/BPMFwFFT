"""
define functions to plot some commonly shown figures
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")

FONTS = {"fontname": "Arial"}

def scater_plot( xd, yd, x_label, y_label, out, 
                figure_size=(3.2, 3.2*6/8), dpi=300, 
                fontsize=8, marker="o", markersize=20, lw=2,
                text_pos=[0.2, 0.8], text_to_title=False):
    """
    """
    reorder = sorted(range(len(xd)), key=lambda i: xd[i])
    xd = [xd[i] for i in reorder]
    yd = [yd[i] for i in reorder]

    coeffs = np.polyfit(xd, yd, 1)
    intercept = coeffs[-1]
    slope = coeffs[-2]
    xl = np.array([np.min(xd), np.max(xd)])
    yl = (slope * xl) + intercept
    line_text = "y = %0.2fx"%slope + (" + " if intercept >= 0 else " - ") + "%0.2f"%np.abs(intercept)

    corr_coef = np.corrcoef([xd, yd])[0][-1]

    rmse = (np.array(xd) - np.array(yd))**2
    rmse = np.sqrt(rmse.mean())

    plt.figure(figsize=figure_size)
    plt.scatter(xd, yd, s=markersize, marker=marker)
    plt.plot(xl, yl, "r", lw=lw)

    ax = plt.axes()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    lower_b_x = np.floor(min(xd))
    upper_b_x = np.ceil(max(xd))

    lower_b_y = np.floor(min(yd))
    upper_b_y = np.ceil(max(yd))

    out_text = "R = %0.2f\nRMSE = %0.2f\n%s"%(corr_coef, rmse, line_text)
    if not text_to_title:
        place_text_at_x = lower_b_x + (upper_b_x - lower_b_x)*text_pos[0]
        place_text_at_y = lower_b_y + (upper_b_y - lower_b_y)*text_pos[1]
        plt.text(place_text_at_x, place_text_at_y, out_text, fontsize=fontsize, **FONTS)
    else:
        out_text = " ".join(out_text.split("\n"))
        plt.title(out_text, fontsize=fontsize, **FONTS)

    plt.xlabel(x_label, fontsize=fontsize, **FONTS)
    plt.ylabel(y_label, fontsize=fontsize, **FONTS)

    axes = plt.gca()
    axes.set_xlim([lower_b_x, upper_b_x])
    axes.set_ylim([lower_b_y, upper_b_y])

    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    return None


def plot_histogram(data, xlabel, out,\
                    title=None, x_range=None, fontsize=8, \
                    figure_size=(3.2, 3.2*6/8), dpi=300):
    plt.figure(figsize=figure_size)
    plt.hist(data)
    if title is not None:
        plt.title(title, fontsize=fontsize, **FONTS)
    plt.xlabel(xlabel, fontsize=fontsize, **FONTS)
    plt.ylabel('Count', fontsize=fontsize, **FONTS)

    ax = plt.axes()
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    if x_range is not None:
        axes = plt.gca()
        axes.set_xlim(x_range)
    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    return None

