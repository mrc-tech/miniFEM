'''
    Questo script prende un modello da JSON e modifica iterativamente le aree degli elementi
    fino a quando non si raggiunge una tensione sigma obiettivo.
'''
import json
import numpy as np
import matplotlib.pyplot as plt

# aggiunge la cartella src al path (cartella dove sta miniFEM)
import sys; import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))

from miniFEM import truss2D as fem


# MODELLO --------------------------------------------------------------------

units = "m, kN, kPa, s"

# importa il modello da un file JSON:
# with open("test_model.json",'r') as f: data = json.load(f)
with open("Long.json",'r') as f: data = json.load(f)
nodes = data['nodes']
elements = data['elements']
constraints = data['constraints']
forces = data['forces']

# --------------------------------------------------------------------

dof = 2 # due spostamenti nel piano (no rotazioni)

def solve_FEM(nodes,elements,constraints,forces):
    K = np.zeros((len(nodes)*dof, len(nodes)*dof)) # matrice di rigidezza del sistema (default: dtype=float64)
    f = fem.create_forces(forces,nodes)
    for e in elements: fem.add_element_stiffness(K, e, nodes) # assembla la matrice di rigidezza K
    for c in constraints: fem.add_constraint(K, c, nodes) # aggiunge i vincoli (sempre alla matrice K)
    u = np.linalg.solve(K,f) # trova gli spostamenti invertendo la matrice K (uguale a "u = np.linalg.inv(K) @ f", ma meno efficiente e meno numericamente stabile)
    
    # salva gli spostamenti nei nodi:
    for i,_ in enumerate(nodes):
        nodes[i]['u'] = u[i*dof+0].item()
        nodes[i]['v'] = u[i*dof+1].item()
    K = np.zeros((len(nodes)*dof, len(nodes)*dof))
    for e in elements: fem.add_element_stiffness(K, e, nodes)
    f = K @ u
    
    # dopo aver calcolato f = K @ u (servono u e v per i nodi)
    element_forces = fem.compute_element_forces(elements, nodes)
    # salva gli sforzi nelle aste:
    N = []
    for i,ef in enumerate(element_forces):
        elements[i]['N'] = ef['N']
        N.append(ef['N'])
    
    return u, f, N

# --------------------------------------------------------------------


sigma_ob = 1  # tensione obiettivo (gli elementi dovrebbero avere tutti questa tensione)
# FARE UNA TENSIONE DIVERSA IN TRAZIONE E IN COMPRESSIONE!!!!
A_min = 0.001 # area minima delle aste
Tol = 0.000000001 # tolleranza per la convergenza
errore = 100*Tol # parametro per uscire dalle iterazioni

count = 0
solve_FEM(nodes,elements,constraints,forces) # risolve per la prima volta
while errore >= Tol: # itera finche' non ha ottimizzato
    for i,_ in enumerate(elements):
        elements[i]['A'] = max([np.abs(elements[i]['N']) / sigma_ob, A_min])
    solve_FEM(nodes,elements,constraints,forces) # risolve con le nuove aree
    count += 1
    errore = 0
    for e in elements:
        errore += pow(np.abs(e['N']/e['A']) - sigma_ob, 2)
    print(f"Iterazione {count}: errore = {errore}")
        

# --------------------------------------------------------------------
# salva il modello:
model = {'model':'Truss 2D', 'units':units, 'author':'Andrea Marchi',  'nodes': nodes, 'elements': elements, 'constraints': constraints, 'forces': forces}
with open('truss2D_model_optimized.json', 'w', encoding='utf-8') as f: json.dump(model, f, ensure_ascii=False, indent='\t')

# PLOT --------------------------------------------------------------------


maxA = max(abs(e['A']) for e in elements) # area massima
minA = min(abs(e['A']) for e in elements) # area minima

fig, ax = plt.subplots()

import matplotlib.patheffects as path_effects
def define_path_effect(**kwargs):
    return [path_effects.Stroke(**kwargs), path_effects.Normal()]

for e in elements:
    i, j = e['i'], e['j']
    A = e['A']

    xi, yi = nodes[i]['x'], nodes[i]['y']
    xj, yj = nodes[j]['x'], nodes[j]['y']
    
    lw = 1 + 6 * abs(A)/maxA

    ax.plot([xi, xj], [yi, yj], color='black', linewidth=lw)

    xm, ym = (xi+xj)/2, (yi+yj)/2
    # ax.text(xm, ym, f"{N:.1f}", color=color, fontsize=9)

    ax.text(
        xm, ym, f"{A:.1f}",
        fontsize=10,
        color='black',
        path_effects=define_path_effect(linewidth=3, foreground="white")
    )

ax.set_aspect('equal')

plt.show()