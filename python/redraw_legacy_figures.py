from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from PIL import Image, ImageDraw, ImageFont
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parent
BASE = Path(__file__).resolve().parents[1] / "data" / "benchmarks"
MEDIA = ROOT / "extracted_media"
OUT = Path(__file__).resolve().parents[1] / "output" / "redraw_legacy_figures"
OUT.mkdir(parents=True, exist_ok=True)
OUT.mkdir(exist_ok=True)

ALGS = ["CPO-DBO", "CPO", "DBO", "IPSO", "GA"]
FILES = {
    "CPO-DBO": BASE / "bench_CPO_DBO_20_paper" / "results_CPO_DBO_20instances.mat",
    "CPO": BASE / "bench_CPO_20_paper" / "results_CPO_20instances.mat",
    "DBO": BASE / "bench_DBO_20_paper" / "results_DBO_20instances.mat",
    "IPSO": BASE / "bench_IPSO_20_paper" / "results_IPSO_20instances.mat",
    "GA": BASE / "bench_GA_20_paper" / "results_GA_20instances.mat",
}
DATA = {name: loadmat(path, squeeze_me=True) for name, path in FILES.items()}

plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"],
                     "mathtext.fontset": "stix", "font.size": 9, "axes.linewidth": 0.8,
                     "axes.grid": False, "figure.facecolor": "white", "axes.facecolor": "white",
                     "pdf.fonttype": 42, "ps.fonttype": 42})
COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7"]


def save(fig, name):
    stem = Path(name).stem
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(OUT / f"{stem}.tiff", dpi=600, bbox_inches="tight", pad_inches=0.03,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


def style_axis(ax):
    ax.grid(False)
    ax.minorticks_off()
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(which="major", direction="in", length=3.5, width=0.8,
                   bottom=True, left=True, top=False, right=False, pad=3)
    ax.tick_params(which="minor", bottom=False, left=False, top=False, right=False)


fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2))
metrics = [
    ("bestDoseRuns", "Best relative exposure cost (a.u.)"),
    ("firstHitIter", "First-hit iteration"),
    ("timeRuns", "Runtime (s)"),
    ("bestLenRuns", "Best route length"),
]
for label, ax, (key, ylabel) in zip("abcd", axes.flat, metrics):
    vals = [np.asarray(DATA[a][key]).ravel() for a in ALGS]
    box = ax.boxplot(vals, tick_labels=ALGS, patch_artist=True, showfliers=False)
    for patch, color in zip(box["boxes"], COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
    for i, values in enumerate(vals, 1):
        ax.scatter(np.full(values.size, i), values, s=10, color=COLORS[i - 1], alpha=0.65)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=20)
    style_axis(ax)
    ax.text(-0.13, 1.02, f"({label})", transform=ax.transAxes, fontsize=11, fontweight="bold")
fig.tight_layout()
save(fig, "Fig8_corrected.png")


fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3))
for alg, color in zip(ALGS, COLORS):
    conv = np.asarray(DATA[alg]["convAll"], float)
    mean = conv.mean(axis=1)
    sd = conv.std(axis=1, ddof=1)
    x = np.arange(1, len(mean) + 1)
    axes[0].plot(x, mean, color=color, lw=1.4, label=alg)
    axes[0].fill_between(x, mean - sd, mean + sd, color=color, alpha=0.10)
    limit = min(1000, len(mean))
    axes[1].plot(x[:limit], mean[:limit], color=color, lw=1.4, label=alg)
    axes[1].fill_between(x[:limit], mean[:limit] - sd[:limit], mean[:limit] + sd[:limit], color=color, alpha=0.10)
axes[0].set_title("(a) Full recorded history")
axes[1].set_title("(b) Early search")
for ax in axes:
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Best relative exposure cost (a.u.)")
    style_axis(ax)
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout()
save(fig, "Fig9_corrected.png")


fig, ax = plt.subplots(figsize=(7.2, 3.5))
for alg, color in zip(ALGS, COLORS):
    values = np.sort(np.asarray(DATA[alg]["bestDoseRuns"]).ravel())
    ax.step(values, np.arange(1, values.size + 1) / values.size, where="post", label=alg, color=color, lw=1.6)
ax.set_xlabel("Best relative exposure cost (a.u.)")
ax.set_ylabel("ECDF")
style_axis(ax)
ax.legend(frameon=False, ncol=2)
fig.tight_layout()
save(fig, "Fig10_corrected.png")


fig, ax = plt.subplots(figsize=(7.2, 3.5))
x = np.arange(len(ALGS))
width = 0.23
for offset, q, label, color in [(-width, 90, "P90", "#0072B2"), (0, 95, "P95", "#D55E00"), (width, 100, "Maximum", "#E69F00")]:
    vals = [np.percentile(np.asarray(DATA[a]["bestDoseRuns"]).ravel(), q) for a in ALGS]
    ax.bar(x + offset, vals, width, label=label, color=color)
ax.set_xticks(x, ALGS)
ax.set_ylabel("Relative exposure cost (a.u.)")
style_axis(ax)
ax.legend(frameon=False)
fig.tight_layout()
save(fig, "Fig11_corrected.png")


FONT_PATH = font_manager.findfont("DejaVu Sans")


def font(size):
    return ImageFont.truetype(FONT_PATH, size=size)


def paste_rotated_label(im, text, xy, size=24, angle=90):
    f = font(size)
    bbox = f.getbbox(text)
    layer = Image.new("RGBA", (bbox[2] - bbox[0] + 12, bbox[3] - bbox[1] + 12), "white")
    draw = ImageDraw.Draw(layer)
    draw.text((6, 6), text, font=f, fill="black")
    layer = layer.rotate(angle, expand=True, fillcolor="white")
    im.paste(layer, xy, layer)


im = Image.open(MEDIA / "image129.jpeg").convert("RGBA")
draw = ImageDraw.Draw(im)
draw.rectangle((1132, 20, 1268, 610), fill="white")
for value, y in [(25, 46), (20, 151), (15, 256), (10, 361), (5, 466), (0, 571)]:
    draw.line((1128, y, 1138, y), fill="black", width=1)
    draw.text((1142, y - 10), str(value), font=font(17), fill="black")
paste_rotated_label(im, "Relative exposure cost (a.u.)", (1195, 190), size=20, angle=90)
im.save(OUT / "Fig4_corrected.png")


im = Image.open(MEDIA / "image131.png").convert("RGBA")
draw = ImageDraw.Draw(im)
for rect in [(55, 60, 105, 300), (715, 60, 765, 300), (55, 415, 105, 670), (715, 415, 765, 670)]:
    draw.rectangle(rect, fill="white")
for xy in [(58, 95), (718, 95), (58, 450), (718, 450)]:
    paste_rotated_label(im, "Relative exposure cost (a.u.)", xy, size=13, angle=90)
for rect, xy in [
    ((530, 23, 650, 45), (550, 25)), ((1190, 23, 1310, 45), (1210, 25)),
    ((530, 378, 650, 400), (550, 380)), ((1190, 378, 1310, 400), (1210, 380)),
]:
    draw.rectangle(rect, fill="white")
    draw.text(xy, "Best cost", font=font(11), fill="black")
im.save(OUT / "Fig6_corrected.png")


im = Image.open(MEDIA / "image140.jpeg").convert("RGBA")
draw = ImageDraw.Draw(im)
draw.rectangle((70, 150, 135, 510), fill="white")
paste_rotated_label(im, "Relative exposure cost (a.u.)", (80, 205), size=17, angle=90)
draw.rectangle((855, 50, 1135, 190), fill="white")
legend_text = ["Run 1: seed 6, best 327.719", "Run 2: seed 7, best 327.719",
               "Run 3: seed 9, best 327.719", "Run 4: seed 15, best 327.719",
               "Run 5: seed 19, best 327.719"]
for i, line in enumerate(legend_text):
    draw.text((880, 60 + 17 * i), line, font=font(12), fill="black")
im.save(OUT / "Fig14_corrected.png")


im = Image.open(MEDIA / "image143.png").convert("RGBA")
draw = ImageDraw.Draw(im)
draw.rectangle((120, 5, 610, 48), fill="white")
draw.rectangle((695, 5, 1250, 48), fill="white")
draw.text((145, 13), "Relative cost = 212.364 a.u., distance = 1046.5", font=font(16), fill="black")
draw.text((715, 13), "Relative cost = 216.782 a.u., distance = 1070.1", font=font(16), fill="black")
im.save(OUT / "Fig16_corrected.png")


im = Image.open(MEDIA / "image144.jpeg").convert("RGBA")
draw = ImageDraw.Draw(im)
draw.rectangle((365, 672, 805, 744), fill="white")
draw.text((425, 690), "Exposure-cost proportion", font=font(24), fill="black")
im.save(OUT / "Fig17_corrected.png")


im = Image.open(MEDIA / "image145.png").convert("RGBA")
draw = ImageDraw.Draw(im)
titles = [
    ((255, 28, 690, 62), (270, 34), "Relative cost = 216.782 a.u., distance = 1070.1"),
    ((875, 28, 1395, 62), (885, 34), "Relative cost = 215.569 a.u., distance = 1060.7"),
    ((255, 588, 690, 622), (270, 594), "Relative cost = 216.777 a.u., distance = 1064.7"),
    ((875, 588, 1395, 622), (885, 594), "Relative cost = 215.411 a.u., distance = 1057.8"),
]
for rect, xy, text in titles:
    draw.rectangle(rect, fill="white")
    draw.text(xy, text, font=font(15), fill="black")
im.save(OUT / "Fig18_corrected.png")

print("corrected figures written to", OUT)
