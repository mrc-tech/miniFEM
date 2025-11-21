#!/usr/bin/env python3
"""
make_truss_tikz.py

Genera un file TikZ che contiene:
- nodi indeformati dichiarati come \node (N1) at (...){};
- nodi deformati dichiarati come \node (D1) at (...){};
- aste indeformate (linee nere)
- aste deformate colorate (rosso = trazione, blu = compressione)
- spessori proporzionali alla forza assiale
"""

import json
import math

# ---------------- CONFIG ----------------
JSON_PATH = "truss2D_model.json"
OUT_TEX = "truss_tikz.tex"

geom_scale = 3.0     # scala della geometria MI DA IDEA CHE QUESTA è SEMPRE 1!!!!!!!!
scale_disp = None    # se None => calcolo automatico
min_linewidth = 0.6
max_linewidth = 3.0
label_forces = True
# ----------------------------------------


def compute_element_forces(data):
    nodes = {n["id"]: n for n in data["nodes"]}
    elements = data["elements"]
    results = []

    # compute max displacement and geometry extent
    xs = [n["x"] for n in data["nodes"]]
    ys = [n["y"] for n in data["nodes"]]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    max_geom_extent = math.hypot(maxx - minx, maxy - miny)
    max_disp = 0.0

    for n in data["nodes"]:
        u = n.get("u", 0.0)
        v = n.get("v", 0.0)
        max_disp = max(max_disp, math.hypot(u, v))

    for el in elements:
        i = el["i"]
        j = el["j"]
        Ni = nodes[i]
        Nj = nodes[j]

        xi, yi = float(Ni["x"]), float(Ni["y"])
        xj, yj = float(Nj["x"]), float(Nj["y"])
        ui, vi = float(Ni.get("u", 0.0)), float(Ni.get("v", 0.0))
        uj, vj = float(Nj.get("u", 0.0)), float(Nj.get("v", 0.0))
        
        L0 = math.hypot(xj-xi, yj-yi)
        extension = ((uj-ui)*(xj-xi) + (vj-vi)*(yj-yi)) / L0
        strain = extension / L0
        
        axial_force = el["E"] * el["A"] * strain

        results.append({
            "i": i, "j": j,
            "xi": xi, "yi": yi, "xj": xj, "yj": yj,
            "ui": ui, "vi": vi, "uj": uj, "vj": vj,
            "L0": L0,
            "axial_force": axial_force
        })

    return results, max_disp, max_geom_extent






# legge il modello:
with open(JSON_PATH, "r", encoding="utf8") as f: data = json.load(f)

# calcola le forze nelle aste:
el_results, max_disp, max_geom_extent = compute_element_forces(data)




# calcola la scala degli spostamenti automaticamente
if scale_disp is None:
    if max_disp <= 0: scale_disp = 1.0
    else: scale_disp = 0.10 * max_geom_extent / max_disp
# calcola la forza massima per scalare lo spessore delle linee
max_force = max(abs(el["axial_force"]) for el in el_results)
if max_force == 0: max_force = 1.0 # evita la divisione per zero
def lw(F): return min_linewidth + (max_linewidth - min_linewidth) * (abs(F) / max_force) # funzione che torna lo spessore della linea in funzione della forza F

out = [] # codice latex (array di stringhe, dove ogni elemento è una riga)

# --- intestazione LaTeX ---
# out.append(r"\documentclass[tikz,border=4pt]{standalone}")
out.append(r"\documentclass[a4paper]{article}")
out.append(r"\usepackage{tikz}")
out.append(r"\begin{document}")
out.append("")
out.append(r"\usetikzlibrary{calc}")
out.append("")

# ---------------- FIGURA 1 (indeformata e deformata) ----------------
out.append("")
out.append(r"\begin{tikzpicture}[scale=1]")

# ---------------- COSTANTI DI VISUALIZZAZIONE ----------------
out.append("% Parametri visualizzazione:")
# out.append(f"\\def\\geomscale{{{geom_scale}}}")
out.append(f"\\def\\dispscale{{{scale_disp}}}")

# ---------------- INDEFORMATA ----------------
out.append("% Nodi indeformati:")
# for n in data["nodes"]: out.append(f"\\node (N{n['id']}) at ({{\\geomscale * {n['x']}}},{{\\geomscale * {n['y']}}}) {{{n['id']}}};")
for n in data["nodes"]: out.append(f"\\node[draw, circle, fill=white, inner sep=2pt, label={{[above right]{n['id']}}}] (N{n['id']}) at ({n['x']},{n['y']}) {{}};")
out.append("% Aste indeformate:")
for el in el_results: out.append(f"\\draw[black, thin] (N{el['i']}) -- (N{el['j']});")


# ---------------- DEFORMATA ----------------
out.append("% Spostamenti nodi:")
for n in data["nodes"]: out.append(f"\\node (u{n['id']}) at ({n['u']},{n['v']}) {{}};")
for n in data["nodes"]: out.append(f"\\node[draw, circle, fill=black, inner sep=2pt] (D{n['id']}) at ($ (N{n['id']}) + \dispscale*(u{n['id']}) $) {{}};")
out.append("% Vettori di spostamento:")
for n in data["nodes"]: out.append(f"\\draw[lightgray,-latex] (N{n['id']}) -- (D{n['id']});")
out.append("% Aste deformate:")
for el in el_results: out.append(f"\\draw[very thick] (D{el['i']}) -- (D{el['j']});")
















out.append(r"\end{tikzpicture}")
out.append("")
# ---------------- FINE FIGURA 1 ----------------


out.append(r"\end{document}")
# ---------------- FINE DOCUMENTO ----------------




# salva il file sul disco:
with open(OUT_TEX, "w", encoding="utf8") as f: f.write("\n".join(out))

# print(f"Creato file TikZ: {OUT_TEX}")


exit()
# ---------------------------------------------------------------------------------------------------------------------------------------------
    

def make_tikz(data, el_results, max_disp, max_geom_extent, out_path,
              geom_scale=1.0, scale_disp=None,
              min_lw=0.6, max_lw=3.0, label_forces=True):

    # compute disp scale automatically
    if scale_disp is None:
        if max_disp <= 0:
            scale_disp = 1.0
        else:
            scale_disp = 0.10 * max_geom_extent / max_disp

    # max force for scaling linewidths
    max_force = max(abs(el["axial_force"]) for el in el_results)
    if max_force == 0:
        max_force = 1.0

    

    # --- begin writing tex ---
    out = []
    out.append(r"\documentclass[tikz,border=4pt]{standalone}")
    out.append(r"\usepackage{tikz}")
    out.append(r"\begin{document}")
    out.append(r"\begin{tikzpicture}[x=1cm,y=1cm]")

    # ---------------- NODI INDEFORMATI ----------------
    out.append("% Nodi indeformati")
    for n in data["nodes"]:
        X = n["x"] * geom_scale
        Y = n["y"] * geom_scale
        nid = n["id"]
        out.append(f"\\node (N{nid}) at ({X:.4f},{Y:.4f}) {{}};")

    # aste indeformate
    out.append("% Aste indeformate (linee nere)")
    for el in el_results:
        i = el["i"]
        j = el["j"]
        out.append(f"\\draw[black, thin] (N{i}) -- (N{j});")

    # ---------------- NODI DEFORMATI ----------------
    out.append("% Nodi deformati")
    for n in data["nodes"]:
        X = (n["x"] + scale_disp * n.get("u", 0)) * geom_scale
        Y = (n["y"] + scale_disp * n.get("v", 0)) * geom_scale
        nid = n["id"]
        out.append(f"\\node (D{nid}) at ({X:.4f},{Y:.4f}) {{}};")

    # ---------------- ASTE DEFORMATE ----------------
    out.append("% Aste deformate (rosso = trazione, blu = compressione)")
    for el in el_results:
        i = el["i"]
        j = el["j"]
        F = el["axial_force"]

        color = "red" if F > 0 else "blue"
        width = lw(F)

        out.append(
            f"\\draw[{color}, line width={width:.3f}pt] (D{i}) -- (D{j});"
        )

        if label_forces:
            mid_i_x = (el["xi"] + scale_disp * el["ui"]) * geom_scale
            mid_i_y = (el["yi"] + scale_disp * el["vi"]) * geom_scale
            mid_j_x = (el["xj"] + scale_disp * el["uj"]) * geom_scale
            mid_j_y = (el["yj"] + scale_disp * el["vj"]) * geom_scale
            xm = 0.5 * (mid_i_x + mid_j_x)
            ym = 0.5 * (mid_i_y + mid_j_y)
            out.append(
                f"\\node[font=\\footnotesize,{color}] at ({xm:.4f},{ym:.4f}) "
                f"{{{F/1000:.3g} kN}};"
            )

    # ---------------- END ----------------
    out.append(r"\end{tikzpicture}")
    out.append(r"\end{document}")

    with open(out_path, "w", encoding="utf8") as f:
        f.write("\n".join(out))

    print(f"Creato file TikZ: {out_path}")



make_tikz(
        data, el_results, max_disp, max_geom_extent,
        OUT_TEX,
        geom_scale=geom_scale,
        scale_disp=scale_disp,
        min_lw=min_linewidth,
        max_lw=max_linewidth,
        label_forces=label_forces
    )