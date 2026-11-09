"""Fig. 2 of GP_manuscript.tex: the workflow diagram of the paper (originally
drawn in Google Drawings, images/Workflow_Figure.png), redrawn as SVG with a
Gaussian process block among the surrogate models.

    python make_workflow_svg.py            # figures/workflow_gp.svg
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
W, H = 1920, 1080
FONT = "font-family:'Times New Roman',Times,serif"
out = []


def add(s):
    out.append(s)


def rect(x, y, w, h, fill, stroke="#595959", sw=2):
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


def text(x, y, s, size=36, anchor="middle", weight="normal", style=""):
    add(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" '
        f'style="{FONT};{style}">{s}</text>')


def lines(x, y, rows, size=36, gap=46):
    for i, r in enumerate(rows):
        text(x, y + i * gap, r, size)


def arrow(x0, y0, x1, y1, w=60, label=None, size=34):
    """Block arrow from (x0, y0) to (x1, y1), horizontal only."""
    d = 1 if x1 > x0 else -1
    head = 45 * d
    body = w / 2
    pts = [(x0, y0 - body), (x1 - head, y0 - body), (x1 - head, y0 - body - 20), (x1, y0),
           (x1 - head, y0 + body + 20), (x1 - head, y0 + body), (x0, y0 + body)]
    add('<polygon points="' + " ".join(f"{a},{b}" for a, b in pts) + '" fill="#efefef" stroke="#595959" stroke-width="2"/>')
    if label:
        text((x0 + x1 - head) / 2, y0 + 12, label, size)


def computer(x, y, s=1.0):
    """Tower PC icon with its top-left corner at (x, y), about 115 x 175 at s = 1."""
    g = f'<g transform="translate({x},{y}) scale({s})">'
    g += '<polygon points="0,25 25,0 115,0 90,25" fill="#7f7f7f" stroke="#333" stroke-width="2"/>'
    g += '<polygon points="90,25 115,0 115,150 90,175" fill="#5f6368" stroke="#333" stroke-width="2"/>'
    g += '<rect x="0" y="25" width="90" height="150" fill="#9aa0a6" stroke="#333" stroke-width="2"/>'
    for yy in (40, 58):
        g += f'<rect x="10" y="{yy}" width="70" height="12" fill="#5f6368" stroke="#333" stroke-width="1"/>'
    g += '<rect x="10" y="95" width="70" height="6" fill="#5f6368"/>'
    g += '<circle cx="45" cy="118" r="5" fill="#7ed957" stroke="#333"/>'
    g += '<circle cx="45" cy="132" r="3" fill="#5f6368"/>'
    g += '<rect x="10" y="145" width="70" height="22" rx="4" fill="#5f6368" stroke="#333" stroke-width="1"/>'
    add(g + "</g>")


add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">')
add(f'<rect width="{W}" height="{H}" fill="white"/>')

# Training data generation
rect(15, 15, 1022, 430, "#c9daf8")
text(35, 82, "Training Data Generation", 46, "start", "bold")
rect(30, 108, 330, 210, "#efefef")
computer(60, 125)
rect(208, 117, 140, 192, "#efefef")
lines(278, 182, ["Single", "edge", "NTM"])
arrow(377, 213, 562, 213, label="Simulate")
rect(573, 92, 375, 242, "#efefef")
text(760, 158, "Single edge datasets", 38)
cols = ["γ", "λ", "ε", "δ", "N<tspan font-size='20' dy='6'>1</tspan>", "N<tspan font-size='20' dy='6'>2</tspan>",
        "φ<tspan font-size='20' dy='6'>τ</tspan>", "Δ<tspan font-size='20' dy='6'>τ</tspan>"]
x0, cw = 588, 43.5
for i, c in enumerate(cols):
    rect(x0 + i * cw, 203, cw, 58, "#efefef", "#999", 1)
    rect(x0 + i * cw, 261, cw, 58, "#efefef", "#999", 1)
    text(x0 + i * cw + cw / 2, 243, c, 26, style="font-style:italic")

# Machine learning surrogate models
rect(1055, 15, 850, 1020, "#d9d2e9")
text(1075, 82, "Machine Learning Surrogate Models", 46, "start", "bold")
px, pw = 1070, 820
# Linear regression
rect(px, 125, pw, 80, "#b4a7d6")
text(px + pw / 2, 180, "Linear regression", 42)
# Symbolic regression
rect(px, 218, pw, 250, "#b4a7d6")
text(px + pw / 2, 268, "Symbolic regression", 42)
nodes = {"+": (1205, 315), "×": (1153, 380), "exp": (1258, 380), "3": (1120, 448), "x1": (1188, 448), "x2": (1258, 448)}
for a, b in [("+", "×"), ("+", "exp"), ("×", "3"), ("×", "x1"), ("exp", "x2")]:
    (xa, ya), (xb, yb) = nodes[a], nodes[b]
    add(f'<line x1="{xa}" y1="{ya + 24}" x2="{xb}" y2="{yb - 24}" stroke="#222" stroke-width="2"/>')
for k, (x, y) in nodes.items():
    rect(x - 30, y - 24, 60, 48, "#efefef", "#595959", 1.5)
    lab = {"x1": "x<tspan font-size='18' dy='5'>1</tspan>", "x2": "x<tspan font-size='18' dy='5'>2</tspan>"}.get(k, k)
    text(x, y + 10, lab, 30 if k != "exp" else 24, style="font-style:italic" if k.startswith("x") else "")
arrow(1305, 380, 1410, 380, w=50)
rect(1425, 320, 445, 120, "#b4a7d6")
text(1647, 395, "3x<tspan font-size='28' dy='8'>1</tspan><tspan dy='-8'> + exp(x</tspan><tspan font-size='28' dy='8'>2</tspan><tspan dy='-8'>)</tspan>", 50)
# Gaussian process
rect(px, 481, pw, 175, "#b4a7d6")
text(px + pw / 2, 528, "Gaussian process", 42)
gx, gy, gw, gh = 1240, 545, 480, 100
add(f'<line x1="{gx}" y1="{gy + gh}" x2="{gx + gw}" y2="{gy + gh}" stroke="#222" stroke-width="2"/>')
add(f'<line x1="{gx}" y1="{gy + gh}" x2="{gx}" y2="{gy}" stroke="#222" stroke-width="2"/>')
import math
xs = [gx + 10 + i * (gw - 20) / 60 for i in range(61)]
mu = [gy + 55 - 28 * math.sin((x - gx) / 75) for x in xs]
obs = [gx + 60, gx + 150, gx + 220, gx + 330, gx + 420]
sd = [4 + 18 * min(1, min(abs(x - o) for o in obs) / 45) for x in xs]
band = [(x, m - s) for x, m, s in zip(xs, mu, sd)] + [(x, m + s) for x, m, s in reversed(list(zip(xs, mu, sd)))]
add('<polygon points="' + " ".join(f"{a:.1f},{b:.1f}" for a, b in band) + '" fill="#6a5acd" fill-opacity="0.35" stroke="none"/>')
add('<polyline points="' + " ".join(f"{a:.1f},{b:.1f}" for a, b in zip(xs, mu)) + '" fill="none" stroke="#3c2a78" stroke-width="3"/>')
for o in obs:
    add(f'<circle cx="{o}" cy="{gy + 55 - 28 * math.sin((o - gx) / 75):.1f}" r="6" fill="#222"/>')
text(gx + gw + 15, gy + 40, "mean", 26, "start")
text(gx + gw + 15, gy + 75, "± s.d.", 26, "start")
# MLP neural network
rect(px, 669, pw, 353, "#b4a7d6")
text(px + pw / 2, 716, "MLP neural network", 42)
text(1480, 760, "Hidden Layers", 32)
add('<polyline points="1320,778 1320,768 1640,768 1640,778" fill="none" stroke="#222" stroke-width="1.5"/>')
text(1243, 800, "Inputs", 32)
inputs = ["γ", "λ", "δ", "ϵ", "N<tspan font-size='12' dy='4'>1</tspan>", "N<tspan font-size='12' dy='4'>2</tspan>"]
iy = [828 + i * 34 for i in range(6)]
hx = [1375, 1468, 1588]
hy = [800 + i * 38 for i in range(4)] + [1000]
out_xy = (1688, 900)
for y in iy:
    for yy in hy:
        add(f'<line x1="1270" y1="{y}" x2="{hx[0] - 16}" y2="{yy}" stroke="#444" stroke-width="0.7"/>')
for yy in hy:
    for yy2 in hy:
        add(f'<line x1="{hx[0] + 16}" y1="{yy}" x2="{hx[1] - 16}" y2="{yy2}" stroke="#444" stroke-width="0.7"/>')
    add(f'<line x1="{hx[2] + 16}" y1="{yy}" x2="{out_xy[0] - 16}" y2="{out_xy[1]}" stroke="#444" stroke-width="0.7"/>')
add(f'<line x1="{out_xy[0] + 16}" y1="{out_xy[1]}" x2="1745" y2="{out_xy[1]}" stroke="#444" stroke-width="0.7"/>')
for i, (y, lab) in enumerate(zip(iy, inputs)):
    add(f'<circle cx="1256" cy="{y}" r="15" fill="#efefef" stroke="#888" stroke-width="1"/>')
    text(1256, y + 7, lab, 20)
for x in hx:
    for yy in hy:
        add(f'<circle cx="{x}" cy="{yy}" r="16" fill="#efefef" stroke="#888" stroke-width="1"/>')
    text(x, 965, "⋮", 30)
text(1528, 905, "…", 30)
add(f'<circle cx="{out_xy[0]}" cy="{out_xy[1]}" r="16" fill="#efefef" stroke="#888" stroke-width="1"/>')

# "Training" cross-shaped connector
add('<polygon points="1020,140 1072,140 1072,195 1130,195 1130,248 1072,248 1072,302 1020,302 1020,248 960,248 960,195 1020,195" '
    'fill="#efefef" stroke="#595959" stroke-width="2"/>')
text(1045, 235, "Training", 36)


# Testing
rect(481, 462, 574, 572, "#d9ead3")
text(502, 527, "Testing", 46, "start", "bold")
rect(509, 561, 466, 212, "#efefef")
computer(552, 580, 0.95)
rect(696, 596, 262, 142, "#efefef")
lines(827, 636, ["Network Wide", "NTM", "simulation"])
rect(509, 784, 466, 212, "#efefef")
computer(552, 805, 0.95)
rect(696, 802, 262, 176, "#efefef")
lines(827, 860, ["Predictive", "accuracy for"])
text(827, 948, "φ<tspan font-size='24' dy='8'>τ</tspan><tspan dy='-8'>, Δ</tspan><tspan font-size='24' dy='8'>τ</tspan>", 40, style="font-style:italic")
arrow(1093, 758, 982, 758)

# Biophysical analysis
rect(15, 462, 443, 572, "#fce5cd")
text(28, 527, "Biophysical Analysis", 46, "start", "bold")
rect(86, 561, 308, 160, "#efefef")
lines(240, 610, ["Anterograde vs.", "retrograde", "transport regimes"])
rect(86, 800, 308, 160, "#efefef")
lines(240, 848, ["Relate Tau flux", "(φ<tspan font-size='24' dy='8'>τ</tspan><tspan dy='-8'>) to NTM</tspan>", "parameters"])
arrow(495, 758, 382, 758)

add("</svg>")
(HERE / "figures").mkdir(exist_ok=True)
(HERE / "figures" / "workflow_gp.svg").write_text("\n".join(out), encoding="utf8")
print("figures/workflow_gp.svg")
