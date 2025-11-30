import numpy as np
import json  # per importare ed esportare i dati del modello
import matplotlib.pyplot as plt  # per plottare il modello


# Le variabili globali vengono viste solo dalle funzioni del modulo
dof = 2 # gradi di liberta' per ogni nodo


# aggiunge la matrice di rigidezza dell'elemento alla matrice di rigidezza K (globale)
# def add_element_stiffness(K: np.ndarray, element: dict[str,float], nodes: dict[str,float]) -> None: # funziona per Python >= 3.9
def add_element_stiffness(K: np.ndarray, element, nodes) -> None:
    i = element['i']; j = element['j'] # id dei nodi di estremità (i, j)
    L = np.sqrt(pow(nodes[j]['x']-nodes[i]['x'],2) + pow(nodes[j]['y']-nodes[i]['y'],2)) # lunghezza dell'elemento
    theta = np.arctan2(nodes[j]['y']-nodes[i]['y'], nodes[j]['x']-nodes[i]['x']) # in gradi: theta_deg = (np.degrees(theta) + 360) % 360
    EA_L = element['E'] * element['A'] / L
    c = np.cos(theta); s = np.sin(theta)
    cc = EA_L*c*c; cs = EA_L*c*s; ss = EA_L*s*s # premoltiplica per risparmiare moltiplicazioni dopo
    i0 = i*dof+0; i1 = i*dof+1; j0 = j*dof+0; j1 = j*dof+1 # i0 rappresenta la componente x del nodo i, mentre i1 la componente y del nodo i
    K[i0,i0] += cc; K[i0,i1] += cs; K[i0,j0] -= cc; K[i0,j1] -= cs
    K[i1,i0] += cs; K[i1,i1] += ss; K[i1,j0] -= cs; K[i1,j1] -= ss
    K[j0,i0] -= cc; K[j0,i1] -= cs; K[j0,j0] += cc; K[j0,j1] += cs
    K[j1,i0] -= cs; K[j1,i1] -= ss; K[j1,j0] += cs; K[j1,j1] += ss

# aggiunge il vincolo
def add_constraint(K: np.ndarray, constraint, nodes) -> None:
    id = constraint['id']
    if constraint['u'] == True: # il grado di libertà orizzontale è vincolato
        for i in range(len(nodes)*dof): K[id*dof+0,i] = 0; K[i,id*dof+0] = 0 # cancella tutte le righe e colonne associate al valore [id*dof+0,id*dof+0]
        K[id*dof+0,id*dof+0] = 1
    if constraint['v'] == True: # il grado di libertà verticale è vincolato
        for i in range(len(nodes)*dof): K[id*dof+1,i] = 0; K[i,id*dof+1] = 0 # cancella tutte le righe e colonne associate al valore [id*dof+1,id*dof+1]
        K[id*dof+1,id*dof+1] = 1
    # DOVREBBE CONTROLLARE CHE ANCHE I RELATIVI VALORI DEL VETTORE f SIANO NULLI!!!!!

# crea il vettore delle forze applicate
# def create_forces(forces, nodes: dict[str,float]) -> np.ndarray: # funziona per Python >= 3.9
def create_forces(forces, nodes) -> np.ndarray:
    f = np.zeros((len(nodes)*dof,1))
    for x in forces: f[x['id']*dof+0], f[x['id']*dof+1] = x['f']
    return f

# calcola le forze interne alle aste
def compute_element_forces(elements, nodes) -> dict:
    forces = []
    for e in elements:
        i = e['i']; j = e['j']

        xi, yi = nodes[i]['x'], nodes[i]['y']
        xj, yj = nodes[j]['x'], nodes[j]['y']

        ui, vi = nodes[i]['u'], nodes[i]['v']
        uj, vj = nodes[j]['u'], nodes[j]['v']

        L = np.hypot(xj - xi, yj - yi)
        c = (xj - xi) / L
        s = (yj - yi) / L

        
        eps = ( (uj - ui)*c + (vj - vi)*s ) / L # deformazione assiale (small displacement)
        sigma = e['E'] * eps # sforzo sigma = E * eps
        N = sigma * e['A'] # forza assiale: N = sigma * A (positiva = trazione)

        forces.append({'element':e, 'N':N, 'sigma':sigma})

    return forces


def plot_element_forces(nodes, element_forces) -> None:
    maxN = max(abs(ef['N']) for ef in element_forces) or 1

    fig, ax = plt.subplots()
    
    import matplotlib.patheffects as path_effects
    def define_path_effect(**kwargs):
        return [path_effects.Stroke(**kwargs), path_effects.Normal()]
    
    for ef in element_forces:
        e = ef['element']
        i, j = e['i'], e['j']
        N = ef['N']

        xi, yi = nodes[i]['x'], nodes[i]['y']
        xj, yj = nodes[j]['x'], nodes[j]['y']

        color = 'blue' if N > 0 else 'red'   # blu trazione, rosso compressione
        lw = 1 + 6 * abs(N)/maxN

        ax.plot([xi, xj], [yi, yj], color=color, linewidth=lw)

        xm, ym = (xi+xj)/2, (yi+yj)/2
        # ax.text(xm, ym, f"{N:.1f}", color=color, fontsize=9)

        ax.text(
            xm, ym, f"{N:.1f}",
            fontsize=10,
            color=color,
            path_effects=define_path_effect(linewidth=3, foreground="white")
        )

    ax.set_aspect('equal')


def plot_forces(nodes, forces) -> None:
    for f in forces:
        n = f['id']
        Fx, Fy = f['f']

        x = nodes[f['id']]['x']
        y = nodes[f['id']]['y']

        scale = 0.2  # scala visualizzazione

        plt.arrow( x-Fx*scale, y-Fy*scale, Fx*scale, Fy*scale,
                  head_width=0.02, head_length=0.04,
                  color='red', length_includes_head=True)
        plt.text(x - Fx*scale/2, y - Fy*scale/2, f"({Fx},{Fy})", color='red')