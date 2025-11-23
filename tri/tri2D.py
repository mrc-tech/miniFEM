import json # per importare ed esportare i dati del modello
import numpy as np
import matplotlib.pyplot as plt # PER IL DEBUG???


# MODELLO --------------------------------------------------------------------

units = "m, kN, kPa, Mg, s"

nodes = [
    {'x':0, 'y':0},
    {'x':1, 'y':0},
    {'x':1, 'y':1}]

elements = [ # definizione dei nodi in senso ANTIORARIO
    {'i':0, 'j':1, 'k':2, 'E':210e3, 'ni':0.3, 't':0.1}]

constraints = [
    {'id':0, 'u':True, 'v':True}, # 'u':True vuol dire che lo spostamento lungo x è vincolato (come OpenSees)
    {'id':1, 'u':True, 'v':True}]

forces = [
    {'id':2, 'f':[1,0]}]

for i,_ in enumerate(nodes): nodes[i]['id'] = i # aggiunge gli 'id' ai nodi (serve solo per facilitare agli umani la lettura degli elementi)

# importa il modello da un file JSON:
# with open("tri2D_model.json",'r') as f: data = json.load(f)
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


# def add_element_stiffness(element): # aggiunge la matrice di rigidezza dell'elemento alla matrice di rigidezza K (globale)
    # i = element['i']; j = element['j']; k = element['k'] # id dei nodi di estremità (i, j, k)
    # xi, yi = nodes[i]['x'], nodes[i]['y']
    # xj, yj = nodes[j]['x'], nodes[j]['y']
    # xk, yk = nodes[k]['x'], nodes[k]['y']
    
    # E = element['E']
    # ni = element['ni']
    # t = element['t'] # spessore dell'elemento
    
    # C = np.array([[1, xi, yi], [1, xj, yj], [1, xk, yk]])
    # # --- CORREZIONE ORIENTAMENTO (la cosa che manca!) ---
    # if np.linalg.det(C) < 0:
        # # swap j e k
        # j, k = k, j
        # xj, yj, xk, yk = xk, yk, xj, yj
        # C = np.array([[1, xi, yi], [1, xj, yj], [1, xk, yk]])
    # IC = np.linalg.inv(C) # inversa della matrice C
    # A = 0.5 * np.linalg.det(C) # area dell'elemento
    # B = np.zeros((3,6))
    # for i in range(3):
        # B[0, 2*i]     = IC[1, i]
        # B[1, 2*i + 1] = IC[2, i]
        # B[2, 2*i]     = IC[2, i]
        # B[2, 2*i + 1] = IC[1, i]
    
    # D = (E/(1-ni*ni)) * np.array([[1, ni, 0], [ni, 1, 0], [0, 0, (1-ni)/2]]) # plane stress
    
    # ke = (B.T @ D @ B) * A*t # matrice di rigidezza dell'elemento (ke = B^T D B det(C)*t/2) ("A @ B" è un modo semplice per "np.matmul(A,B)")
    
    # i0 = i*dof+0; i1 = i*dof+1; j0 = j*dof+0; j1 = j*dof+1 # i0 rappresenta la componente x del nodo i, mentre i1 la componente y del nodo i
    # K[i0,i0] += ke[0,0]; K[i0,i1] += ke[0,1]; K[i0,j0] += ke[0,2]; K[i0,j1] += ke[0,3]
    # K[i1,i0] += ke[1,0]; K[i1,i1] += ke[1,1]; K[i1,j0] += ke[1,2]; K[i1,j1] += ke[1,3]
    # K[j0,i0] += ke[2,0]; K[j0,i1] += ke[2,1]; K[j0,j0] += ke[2,2]; K[j0,j1] += ke[2,3]
    # K[j1,i0] += ke[3,0]; K[j1,i1] += ke[3,1]; K[j1,j0] += ke[3,2]; K[j1,j1] += ke[3,3]

def add_element_stiffness_2(element):
    """
    Aggiunge la matrice di rigidezza dell'elemento triangolare 2D (CST) alla matrice globale K
    """
    i, j, k = element['i'], element['j'], element['k'] # nodi dell'elemento
    xi, yi = nodes[i]['x'], nodes[i]['y']
    xj, yj = nodes[j]['x'], nodes[j]['y']
    xk, yk = nodes[k]['x'], nodes[k]['y']
    
    E = element['E']
    ni = element['ni']
    t = element['t'] # spessore dell'elemento
    
    A = abs((xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj)) / 2) # area del triangolo (sempre positiva)
    if A == 0: raise ValueError("Triangolo degenerato con area zero!")
    
    # --- coefficienti b e c ---
    b1, b2, b3 = yj - yk, yk - yi, yi - yj
    c1, c2, c3 = xk - xj, xi - xk, xj - xi
    B = (1/(2*A)) * np.array([
        [b1,  0, b2,  0, b3,  0],
        [ 0, c1,  0, c2,  0, c3],
        [c1, b1, c2, b2, c3, b3]
    ])
    
    D = (E/(1-ni*ni)) * np.array([[1, ni, 0], [ni, 1, 0], [0, 0, (1-ni)/2]]) # plane stress
    
    ke = B.T @ D @ B * t * A # matrice di rigidezza dell'elemento (ke = B^T D B * A*t) ("A @ B" è un modo semplice per "np.matmul(A,B)")
    
    # assemblaggio nella matrice globale K
    dof_map = [i*dof, i*dof+1, j*dof, j*dof+1, k*dof, k*dof+1]
    for r in range(6):
        for c in range(6):
            K[dof_map[r], dof_map[c]] += ke[r, c]

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
for e in elements: add_element_stiffness_2(e) # assembla la matrice di rigidezza K
for c in constraints: add_constraint(c) # aggiunge i vincoli (sempre alla matrice K)
u = np.linalg.solve(K,f) # trova gli spostamenti invertendo la matrice K (uguale a "u = np.linalg.inv(K) @ f", ma meno efficiente e meno numericamente stabile)


# salva gli spostamenti nei nodi:
for i,_ in enumerate(nodes):
    nodes[i]['u'] = u[i*dof+0].item()
    nodes[i]['v'] = u[i*dof+1].item()


print('spostamenti:\n', u)
# ASSEMBLA DI NUOVO LA MATRICE K PER CALCOLARE LE FORZE....
K = np.zeros((N*dof, N*dof))
for e in elements: add_element_stiffness_2(e)
forze = K @ u
print('forze:\n', forze)

# CALCOLARE LE FORZE NELLE ASTE........




# exit()
# --------------------------------------------------------------------
# salva il modello:
model = {'model':'Tri 2D', 'units':units, 'author':'Andrea Marchi',  'nodes': nodes, 'elements': elements, 'constraints': constraints, 'forces': forces}
with open('tri2D_model.json', 'w', encoding='utf-8') as f: json.dump(model, f, ensure_ascii=False, indent='\t')

# PLOT -------------------------------------------------------------------- per il debug??..........................

# indeformata:
for e in elements: plt.plot([nodes[e['i']]['x'], nodes[e['j']]['x'], nodes[e['k']]['x'], nodes[e['i']]['x']], [nodes[e['i']]['y'], nodes[e['j']]['y'], nodes[e['k']]['y'], nodes[e['i']]['y']], '-ok')
# STAMPARE ANCHE GLI ASSI LOCALI (almeno il vettore che va da "i" a "j")
# deformata:
if 'u' in nodes[0]: # controlla se ci sono gli spostamenti da mostrare
    mul = 500 # scala che moltiplica gli spostamenti (multiplier)
    for e in elements:
        plt.plot([nodes[e['i']]['x']+mul*nodes[e['i']]['u'], nodes[e['j']]['x']+mul*nodes[e['j']]['u'], nodes[e['k']]['x']+mul*nodes[e['k']]['u'], nodes[e['i']]['x']+mul*nodes[e['i']]['u']], [nodes[e['i']]['y']+mul*nodes[e['i']]['v'], nodes[e['j']]['y']+mul*nodes[e['j']]['v'], nodes[e['k']]['y']+mul*nodes[e['k']]['v'], nodes[e['i']]['y']+mul*nodes[e['i']]['v']], '-or')
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



