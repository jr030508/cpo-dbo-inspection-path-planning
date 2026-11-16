import csv
import heapq
import json
import math
import os
import sys
import time
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[key] = "1"

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from matplotlib.patches import Rectangle

BASE = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE))
from solver_suite import run_cpo_dbo  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output" / Path(__file__).stem
OUT.mkdir(parents=True, exist_ok=True)
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)

NODES = np.array([
    [120,0],[131,18],[154,11],[164,24],[174,35],[193,48],[187,67],[165,73],[110,16],
    [73,12],[55,29],[72,46],[33,52],[13,64],[29,81],[62,79],[175,91],[198,91],
    [189,99],[173,106],[53,92],[69,104],[98,107],[61,119],[87,127],[186,118],
    [172,151],[160,139],[156,150],[143,123],[131,134],[43,135],[12,143],[52,152],
    [38,161],[22,172],[40,192],[123,143],[135,166],[140,172],[140,183],[109,173],
    [96,193],[87,179],[73,166],[64,188],[83,25],[42,52],[18,88],[23,156],[117,181]
], dtype=float)
SOURCES = np.array([[35,35],[30,110],[88,70],[88,142],[138,107],[160,50],[168,180]], float)
STRENGTH = np.array([8000,12000,10000,15000,9000,13000,11000], float)
OBSTACLES = np.array([[55,12,14,10],[96,22,18,14],[5,115,14,14],[110,150,16,18],[150,120,18,14]], float)
MARGIN = 2.0
GRID = 2.0
DS = 1.0
EPS_LEN = 1e-5


def radiation(x, y):
    value = np.zeros(np.broadcast(x, y).shape, float)
    for (sx, sy), q in zip(SOURCES, STRENGTH):
        value += q / ((x-sx)**2 + (y-sy)**2 + 10.0)
    return value


def line_integral(p1, p2):
    length = float(np.linalg.norm(p2-p1))
    n = max(2, int(math.ceil(length/DS)) + 1)
    t = np.linspace(0.0, 1.0, n)
    x = p1[0] + t*(p2[0]-p1[0])
    y = p1[1] + t*(p2[1]-p1[1])
    dose = float(np.trapezoid(radiation(x, y), dx=length/(n-1)) / 100.0)
    return length, dose


def segment_hits_rect(p1, p2, rect):
    xmin, ymin, width, height = rect
    xmin -= MARGIN; ymin -= MARGIN
    xmax = xmin + width + 2*MARGIN; ymax = ymin + height + 2*MARGIN
    dx, dy = p2-p1
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, p1[0]-xmin), (dx, xmax-p1[0]),
                 (-dy, p1[1]-ymin), (dy, ymax-p1[1])):
        if abs(p) < 1e-15:
            if q < 0:
                return False
        else:
            r = q/p
            if p < 0:
                if r > t1: return False
                if r > t0: t0 = r
            else:
                if r < t0: return False
                if r < t1: t1 = r
    return t0 <= t1


def has_los(p1, p2, obstacles):
    return not any(segment_hits_rect(p1, p2, rect) for rect in obstacles)


def build_grid(obstacles):
    xvec = np.arange(0.0, 200.0 + GRID/2, GRID)
    yvec = np.arange(0.0, 200.0 + GRID/2, GRID)
    xx, yy = np.meshgrid(xvec, yvec)
    occ = np.zeros(xx.shape, bool)
    for xmin, ymin, width, height in obstacles:
        occ |= ((xx >= xmin-MARGIN) & (xx <= xmin+width+MARGIN) &
                (yy >= ymin-MARGIN) & (yy <= ymin+height+MARGIN))
    rad = radiation(xx, yy)
    return xvec, yvec, occ, rad


def node_rc(point):
    return int(round(point[1]/GRID)), int(round(point[0]/GRID))


def astar(start, goal, occ, rad):
    sr, sc = node_rc(start); gr, gc = node_rc(goal)
    moves = [(-1,0,1.0),(1,0,1.0),(0,-1,1.0),(0,1,1.0),
             (-1,-1,math.sqrt(2)),(-1,1,math.sqrt(2)),(1,-1,math.sqrt(2)),(1,1,math.sqrt(2))]
    h_per_m = float(rad.min()/100.0 + EPS_LEN)
    g = {(sr, sc): 0.0}
    parent = {}
    heap = [(h_per_m*GRID*math.hypot(gr-sr, gc-sc), 0.0, sr, sc)]
    closed = set()
    while heap:
        _, current_g, r, c = heapq.heappop(heap)
        if (r,c) in closed: continue
        if (r,c) == (gr,gc): break
        closed.add((r,c))
        for dr, dc, mult in moves:
            nr, nc = r+dr, c+dc
            if nr < 0 or nr >= occ.shape[0] or nc < 0 or nc >= occ.shape[1]: continue
            if occ[nr,nc] or (nr,nc) in closed: continue
            step_len = mult*GRID
            step_dose = step_len*rad[nr,nc]/100.0
            ng = current_g + step_dose + EPS_LEN*step_len
            if ng < g.get((nr,nc), float("inf")):
                g[(nr,nc)] = ng
                parent[(nr,nc)] = (r,c)
                h = h_per_m*GRID*math.hypot(gr-nr, gc-nc)
                heapq.heappush(heap, (ng+h, ng, nr, nc))
    if (gr,gc) not in g:
        raise RuntimeError(f"A* failed {start}->{goal}")
    path = [(gr,gc)]
    while path[-1] != (sr,sc):
        path.append(parent[path[-1]])
    path.reverse()
    xy = np.array([[c*GRID, r*GRID] for r,c in path], float)
    length = 0.0; dose = 0.0
    for (r0,c0),(r1,c1) in zip(path[:-1], path[1:]):
        step = GRID*math.hypot(r1-r0, c1-c0)
        length += step
        dose += step*rad[r1,c1]/100.0
    return xy, length, dose


def build_matrices(active_ids):
    obstacles = OBSTACLES[np.array(active_ids, int)-1] if active_ids else np.empty((0,4))
    _, _, occ, rad = build_grid(obstacles)
    n = len(NODES)
    length = np.zeros((n,n)); dose = np.zeros((n,n)); edge_type = np.zeros((n,n), int)
    paths = {}
    for i in range(n):
        for j in range(i+1,n):
            if has_los(NODES[i], NODES[j], obstacles):
                le, de = line_integral(NODES[i], NODES[j])
                path = np.vstack([NODES[i], NODES[j]])
            else:
                path, le, de = astar(NODES[i], NODES[j], occ, rad)
                edge_type[i,j] = edge_type[j,i] = 1
            length[i,j] = length[j,i] = le
            dose[i,j] = dose[j,i] = de
            paths[(i,j)] = path
    return obstacles, rad, length, dose, dose + EPS_LEN*length, edge_type, paths


def solve_case(label, active_ids):
    started = time.perf_counter()
    obstacles, rad, length, dose, cost, edge_type, paths = build_matrices(active_ids)
    build_seconds = time.perf_counter()-started
    result = run_cpo_dbo(cost, 6, m=400, tscale=30000.0, elite_ratio=0.02,
                         twoopt_attempts=1, budget=1_960_400)
    route = np.r_[0, result[1], 0]
    total_dose = float(sum(dose[a,b] for a,b in zip(route[:-1], route[1:])))
    total_length = float(sum(length[a,b] for a,b in zip(route[:-1], route[1:])))
    astar_edges = int(sum(edge_type[a,b] for a,b in zip(route[:-1], route[1:])))
    astar_dose = float(sum(dose[a,b] for a,b in zip(route[:-1], route[1:]) if edge_type[a,b]))
    astar_length = float(sum(length[a,b] for a,b in zip(route[:-1], route[1:]) if edge_type[a,b]))
    return {
        "label": label, "active_obstacles": active_ids, "obstacles": obstacles, "rad": rad,
        "length_matrix": length, "dose_matrix": dose, "cost_matrix": cost,
        "edge_type": edge_type, "paths": paths, "route": route,
        "total_dose": total_dose, "total_length": total_length, "astar_edges": astar_edges,
        "astar_dose_share_pct": 100*astar_dose/total_dose,
        "astar_length_share_pct": 100*astar_length/total_length,
        "matrix_build_seconds": build_seconds,
    }


def style_axis(ax):
    ax.grid(False); ax.minorticks_off()
    for side in ("top","right","bottom","left"):
        ax.spines[side].set_visible(True); ax.spines[side].set_linewidth(0.8)
    ax.tick_params(which="major", direction="in", length=3.5, width=0.8,
                   bottom=True, left=True, top=False, right=False, pad=3)
    ax.tick_params(which="minor", bottom=False, left=False, top=False, right=False)


def draw_case(ax, case):
    xx, yy = np.meshgrid(np.arange(0,201,GRID), np.arange(0,201,GRID))
    ax.contourf(xx, yy, case["rad"], levels=28, cmap="cividis", alpha=0.82)
    for rect in case["obstacles"]:
        ax.add_patch(Rectangle((rect[0],rect[1]),rect[2],rect[3], facecolor="#8C8C8C",
                               edgecolor="black", lw=0.7, zorder=4))
    route = case["route"]
    for a,b in zip(route[:-1], route[1:]):
        key = (min(a,b), max(a,b)); path = case["paths"][key]
        if a > b: path = path[::-1]
        if case["edge_type"][a,b]:
            ax.plot(path[:,0],path[:,1],color="#009E73",lw=1.4,ls="--",zorder=5)
        else:
            ax.plot(path[:,0],path[:,1],color="#D55E00",lw=1.05,zorder=5)
    ax.scatter(NODES[:,0],NODES[:,1],s=10,c="black",edgecolors="white",linewidths=0.25,zorder=7)
    for idx,(x,y) in enumerate(NODES,1):
        ax.text(x+1.2,y+1.2,str(idx),fontsize=5.0,color="white",zorder=8,
                path_effects=[pe.withStroke(linewidth=1.2,foreground="black")])
    ax.scatter(SOURCES[:,0],SOURCES[:,1],marker="*",s=48,c="#E69F00",
               edgecolors="black",linewidths=0.4,zorder=8)
    ax.set_xlim(0,200); ax.set_ylim(0,200); ax.set_aspect("equal")
    ax.set_xticks(np.arange(0,201,50)); ax.set_yticks(np.arange(0,201,50))
    style_axis(ax)


def save_all(fig, name):
    fig.savefig(FIG/f"{name}.png",dpi=600,bbox_inches="tight",pad_inches=0.03)
    fig.savefig(FIG/f"{name}.tiff",dpi=600,bbox_inches="tight",pad_inches=0.03,
                pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(FIG/f"{name}.pdf",bbox_inches="tight",pad_inches=0.03)
    fig.savefig(FIG/f"{name}.svg",bbox_inches="tight",pad_inches=0.03)


def main():
    mpl.rcParams.update({"font.family":"serif","font.serif":["Times New Roman"],
                         "mathtext.fontset":"stix","font.size":8.5,"axes.titlesize":8.5,
                         "axes.labelsize":9,"xtick.labelsize":7,"ytick.labelsize":7,
                         "legend.fontsize":8,"axes.grid":False,"figure.facecolor":"white",
                         "axes.facecolor":"white","pdf.fonttype":42,"ps.fonttype":42})
    configs = [
        ("Obstacle-free", []),
        ("1,2,3,4,5", [1,2,3,4,5]),
        ("2,3,4,5", [2,3,4,5]),
        ("1,3,4,5", [1,3,4,5]),
        ("3,4,5", [3,4,5]),
    ]
    cases = [solve_case(label, ids) for label,ids in configs]
    summary = []
    for c in cases:
        summary.append({k:c[k] for k in ("label","active_obstacles","total_dose","total_length",
                                         "astar_edges","astar_dose_share_pct","astar_length_share_pct",
                                         "matrix_build_seconds")})
        np.savez_compressed(OUT/f"obstacle_case_{c['label'].replace(',','_')}.npz",
                            route=c["route"], dose=c["dose_matrix"], length=c["length_matrix"],
                            cost=c["cost_matrix"], edge_type=c["edge_type"])
    with (OUT/"standardized_obstacle_results.csv").open("w",newline="",encoding="utf-8-sig") as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0])); w.writeheader(); w.writerows(summary)
    (OUT/"standardized_obstacle_results.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")

    # Fig. 20 compares the obstacle-free and full-obstacle tours generated by
    # the same standardized solver and parameter set.
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.55))
    fig.subplots_adjust(left=.08,right=.985,bottom=.14,top=.69,wspace=.15)
    for i,(ax,c) in enumerate(zip(axes,cases[:2])):
        draw_case(ax,c)
        ax.set_title(f"Cost = {c['total_dose']:.2f} a.u.; distance = {c['total_length']:.1f}",pad=6)
        ax.text(-.10,1.10,f"({chr(97+i)})",transform=ax.transAxes,fontweight="bold",fontsize=10,
                ha="left",va="bottom")
    handles=[
        mpl.lines.Line2D([],[],color="black",marker="s",mfc="#8C8C8C",lw=0,ms=5,label="Obstacle"),
        mpl.lines.Line2D([],[],color="black",marker="o",lw=0,ms=4,label="Inspection node"),
        mpl.lines.Line2D([],[],color="#E69F00",marker="*",lw=0,ms=7,label="Radiation source"),
        mpl.lines.Line2D([],[],color="#D55E00",lw=1.4,label="LOS edge"),
        mpl.lines.Line2D([],[],color="#009E73",lw=1.4,ls="--",label="A* detour edge"),
    ]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.99),ncol=3,frameon=False,
               handlelength=2,columnspacing=1.6)
    fig.text(.52,.03,"X coordinate (m)",ha="center",fontsize=9)
    fig.text(.018,.39,"Y coordinate (m)",va="center",rotation=90,fontsize=9)
    save_all(fig,"Fig20_path_comparison_standardized"); plt.close(fig)

    # Fig. 21 from the same standardized full-obstacle run.
    full = cases[1]
    los = np.array([100-full["astar_dose_share_pct"],100-full["astar_length_share_pct"]])
    ast = 100-los
    fig,ax=plt.subplots(figsize=(3.55,3.05)); fig.subplots_adjust(left=.18,right=.97,bottom=.19,top=.70)
    x=np.arange(2)
    ax.bar(x,los,.52,color="#0072B2",edgecolor="black",lw=.6,label="LOS edges")
    ax.bar(x,ast,.52,bottom=los,color="#E69F00",edgecolor="black",lw=.6,label="A* detour edges")
    for i in range(2):
        ax.text(i,los[i]/2,f"{los[i]:.1f}%",ha="center",va="center",color="white")
        ax.text(i,los[i]+ast[i]/2,f"{ast[i]:.1f}%",ha="center",va="center",color="black")
    ax.set_xticks(x,["Exposure cost","Distance"]); ax.set_ylabel("Share of tour total (%)")
    ax.set_ylim(0,105); ax.set_title("Edge composition of the standardized obstacle-aware tour",pad=7)
    style_axis(ax)
    fig.legend(loc="upper center",bbox_to_anchor=(.57,.985),ncol=2,frameon=False)
    save_all(fig,"Fig21_edge_composition_standardized"); plt.close(fig)

    # Fig. 22 fully redrawn from standardized matrices and tours.
    fig,axes=plt.subplots(2,2,figsize=(7.2,6.35))
    fig.subplots_adjust(left=.08,right=.985,bottom=.075,top=.77,wspace=.16,hspace=.30)
    panels=cases[1:]
    for i,(ax,c) in enumerate(zip(axes.ravel(),panels)):
        draw_case(ax,c)
        ax.set_title(f"Cost = {c['total_dose']:.2f} a.u.; distance = {c['total_length']:.1f}; A* edges = {c['astar_edges']}",pad=6)
        ax.text(-.10,1.10,f"({chr(97+i)})",transform=ax.transAxes,fontweight="bold",fontsize=10,
                ha="left",va="bottom")
    handles=[
        mpl.lines.Line2D([],[],color="black",marker="s",mfc="#8C8C8C",lw=0,ms=5,label="Obstacle"),
        mpl.lines.Line2D([],[],color="black",marker="o",lw=0,ms=4,label="Inspection node"),
        mpl.lines.Line2D([],[],color="#E69F00",marker="*",lw=0,ms=7,label="Radiation source"),
        mpl.lines.Line2D([],[],color="#D55E00",lw=1.4,label="LOS edge"),
        mpl.lines.Line2D([],[],color="#009E73",lw=1.4,ls="--",label="A* detour edge"),
    ]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.99),ncol=3,frameon=False,
               handlelength=2,columnspacing=1.6)
    fig.text(.52,.018,"X coordinate (m)",ha="center",fontsize=9)
    fig.text(.018,.42,"Y coordinate (m)",va="center",rotation=90,fontsize=9)
    save_all(fig,"Fig22_obstacle_combinations_standardized"); plt.close(fig)
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == "__main__":
    main()
