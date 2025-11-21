import json # per importare ed esportare i dati del modello
import numpy as np
import matplotlib.pyplot as plt # PER IL DEBUG???


# MODELLO --------------------------------------------------------------------

units = "m, kN, kPa, s"

nodes = [
    {'x':0, 'y':0},
    {'x':1, 'y':0},
    {'x':0, 'y':1}]

elements = [
    {'i':0, 'j':1, 'E':210e3, 'A':0.1},
    {'i':1, 'j':2, 'E':210e3, 'A':0.1},
    {'i':0, 'j':2, 'E':210e3, 'A':0.1}]

constraints = [
    {'id':0, 'u':True, 'v':True}, # 'u':True vuol dire che lo spostamento lungo x è vincolato (come OpenSees)
    {'id':1, 'u':False, 'v':True}]

forces = [
    {'id':2, 'f':[1,0]}]

for i,_ in enumerate(nodes): nodes[i]['id'] = i # aggiunge gli 'id' ai nodi (serve solo per facilitare agli umani la lettura degli elementi)

# importa il modello da un file JSON:
# with open("truss2D_model.json",'r') as f: data = json.load(f)
# nodes = data['nodes']
# elements = data['elements']
# constraints = data['constraints']
# forces = data['forces']

# --------------------------------------------------------------------

N = len(nodes) # numero di nodi
dof = 2 # due spostamenti nel piano (no rotazioni)


# print('nodes:'); print(nodes)
# print('elements:'); print(elements)

# matrice di rigidezza del sistema
K = np.zeros((N*dof, N*dof)) # default: dtype=float64
# FARE CON LE MATRICI SPARSE (se il numero di nodi è eccessivo)


def add_element_stiffness(element): # aggiunge la matrice di rigidezza dell'elemento alla matrice di rigidezza K (globale)
    i = element['i']; j = element['j'] # id dei nodi di estremità (i, j)
    L = np.sqrt(pow(nodes[j]['x']-nodes[i]['x'],2) + pow(nodes[j]['y']-nodes[i]['y'],2)) # lunghezza dell'elemento
    # theta = np.arctan((nodes[j]['y']-nodes[i]['y'])/(nodes[j]['x']-nodes[i]['x'])) # QUESTO DA PROBLEMI IN CASO DI ELEMENTI VERTICALI
    theta = np.arctan2(nodes[j]['y']-nodes[i]['y'], nodes[j]['x']-nodes[i]['x']) # in gradi: theta_deg = (np.degrees(theta) + 360) % 360
    # print(element)
    # print(L)
    # print((np.degrees(theta) + 360) % 360)
    k = element['E'] * element['A'] / L
    c = np.cos(theta); s = np.sin(theta)
    cc = k*c*c; cs = k*c*s; ss = k*s*s # premoltiplica per risparmiare moltiplicazioni dopo
    i0 = i*dof+0; i1 = i*dof+1; j0 = j*dof+0; j1 = j*dof+1 # i0 rappresenta la componente x del nodo i, mentre i1 la componente y del nodo i
    #      i0  i1  j0  j1
    # i0   cc, cs,-cc,-cs
    # i1   cs, ss,-cs,-ss
    # j0  -cc,-cs, cc, cs
    # j1  -cs,-ss, cs, ss
    K[i0,i0] += cc; K[i0,i1] += cs; K[i0,j0] -= cc; K[i0,j1] -= cs
    K[i1,i0] += cs; K[i1,i1] += ss; K[i1,j0] -= cs; K[i1,j1] -= ss
    K[j0,i0] -= cc; K[j0,i1] -= cs; K[j0,j0] += cc; K[j0,j1] += cs
    K[j1,i0] -= cs; K[j1,i1] -= ss; K[j1,j0] += cs; K[j1,j1] += ss

def add_constraint(constraint): # aggiunge il vincolo
    id = constraint['id']
    if constraint['u'] == True: # il grado di libertà orizzontale è vincolato
        for i in range(N*dof): K[id*dof+0,i] = 0; K[i,id*dof+0] = 0 # cancella tutte le righe e colonne associate al valore [id*dof+0,id*dof+0]
        K[id*dof+0,id*dof+0] = 1
    if constraint['v'] == True: # il grado di libertà verticale è vincolato
        for i in range(N*dof): K[id*dof+1,i] = 0; K[i,id*dof+1] = 0 # cancella tutte le righe e colonne associate al valore [id*dof+1,id*dof+1]
        K[id*dof+1,id*dof+1] = 1
    # DOVREBBE CONTROLLARE CHE ANCHE I RELATIVI VALORI DEL VETTORE f SIANO NULLI!!!!!

def create_forces(forces): # crea il vettore delle forze applicate
    f = np.zeros((N*dof,1))
    for x in forces: f[x['id']*dof+0], f[x['id']*dof+1] = x['f']
    return f


f = create_forces(forces)
for e in elements: add_element_stiffness(e) # assembla la matrice di rigidezza K
for c in constraints: add_constraint(c) # aggiunge i vincoli (sempre alla matrice K)
u = np.linalg.solve(K,f) # trova gli spostamenti invertendo la matrice K (uguale a "u = np.linalg.inv(K) @ f", ma meno efficiente e meno numericamente stabile)


# salva gli spostamenti nei nodi:
for i,_ in enumerate(nodes):
    nodes[i]['u'] = u[i*dof+0].item()
    nodes[i]['v'] = u[i*dof+1].item()


print('spostamenti:\n', u)
# ASSEMBLA DI NUOVO LA MATRICE K PER CALCOLARE LE FORZE....
K = np.zeros((N*dof, N*dof))
for e in elements: add_element_stiffness(e)
forze = K @ u
print('forze:\n', forze)

# CALCOLARE LE FORZE NELLE ASTE........




# exit()
# --------------------------------------------------------------------
# salva il modello:
model = {'model':'Truss 2D', 'units':units, 'author':'Andrea Marchi',  'nodes': nodes, 'elements': elements, 'constraints': constraints, 'forces': forces}
with open('truss2D_model.json', 'w', encoding='utf-8') as f: json.dump(model, f, ensure_ascii=False, indent='\t')

# PLOT -------------------------------------------------------------------- per il debug??..........................

# indeformata:
for e in elements: plt.plot([nodes[e['i']]['x'], nodes[e['j']]['x']], [nodes[e['i']]['y'], nodes[e['j']]['y']], '-ok')
# STAMPARE ANCHE GLI ASSI LOCALI (almeno il vettore che va da "i" a "j")
# deformata:
if 'u' in nodes[0]: # controlla se ci sono gli spostamenti da mostrare
    mul = 500 # scala che moltiplica gli spostamenti (multiplier)
    for e in elements:
        plt.plot([nodes[e['i']]['x']+mul*nodes[e['i']]['u'], nodes[e['j']]['x']+mul*nodes[e['j']]['u']],
                 [nodes[e['i']]['y']+mul*nodes[e['i']]['v'], nodes[e['j']]['y']+mul*nodes[e['j']]['v']], '-or')
    for i,n in enumerate(nodes):
        scale = 0.2 # scala visualizzazione
        plt.arrow(n['x']+mul*n['u'], n['y']+mul*n['v'], forze[i*dof+0].item()*scale, forze[i*dof+1].item()*scale, head_width=0.02, head_length=0.04, color='red', length_includes_head=True)
        # plt.text(nodes[f['id']]['x'] + f['f'][0]*scale, nodes[f['id']]['y'] + f['f'][1]*scale, f"({f['f'][0]},{f['f'][1]})", color='red')
plt.show()


# exit()

# # stampa la matrice K
# plt.imshow(K, cmap='viridis', origin='lower')
# for i in range(K.shape[0]):
    # for j in range(K.shape[1]):
        # plt.text(j, i, str(K[i, j]), ha='center', va='center', color='white') # mostra anche i valori
# plt.colorbar()  # aggiunge la scala dei valori
# plt.show()




#........................................................................
# CODICE AGGIUNTO DA CHATGPT..................................



def compute_element_forces(elements, nodes):
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

        # deformazione assiale (small displacement)
        eps = ( (uj - ui)*c + (vj - vi)*s ) / L

        # sforzo sigma = E * eps
        sigma = e['E'] * eps

        # forza assiale: N = sigma * A (positiva = trazione)
        N = sigma * e['A']

        forces.append({'element':e, 'N':N, 'sigma':sigma})

    return forces


# dopo aver calcolato f = K @ u (servono u e v per i nodi)
element_forces = compute_element_forces(elements, nodes)

print("\nFORZE NELLE ASTE:")
for ef in element_forces:
    i = ef['element']['i']; j = ef['element']['j']
    print(f"Elemento {i}-{j}:  N = {ef['N']:.3f}   sigma = {ef['sigma']:.3f}")

def plot_forces(nodes, forces):
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


import matplotlib.patheffects as path_effects
def define_path_effect(**kwargs):
    return [path_effects.Stroke(**kwargs), path_effects.Normal()]



def plot_element_forces(nodes, element_forces):
    maxN = max(abs(ef['N']) for ef in element_forces) or 1

    fig, ax = plt.subplots()
    
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
    


plot_element_forces(nodes, element_forces)
plot_forces(nodes, forces)

# VORREI FARE ANCHE CHE SALVA GLI SFORZI SULLE ASTE NEL JSON...
# ANCHE LE FORZE NEI NODI


plt.show()