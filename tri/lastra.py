'''
    modello di un pannello in CLS soggetto ad una forza verticale e una orizzontale (tutto nel piano).
    Cerco di capire la resistenza del pannello quando vedo che lo sforzo di Von Mises
    raggiunge la resistenza a compressione
    (SUPER APPROSSIMAZIONE VISTO CHE SI CONSIDERA RESISTENTE ANCHE A TRAZIONE)
    
    con una forza verticale di 400 kN si raggiunge una sigma di 50 MPa con F_hor = 1340 kN (Nx=20*2; Ny=30*2)
'''
import json # per importare ed esportare i dati del modello
import numpy as np
from scipy.linalg import eigh # per l'analisi modale
import matplotlib.pyplot as plt # PER IL DEBUG???
import matplotlib.tri as tri


# MODELLO --------------------------------------------------------------------

units = "m, kN, kPa, Mg, s"

B = 2.4 # base
H = 3.1 # altezza
Nx = 20*2 # numero di suddivisioni lungo x (con 20*5 ci mette qualche minuto e occupa più di 10Gb di RAM)
Ny = 30*2 # numero di suddivisioni lungo y
Delta_x = B/Nx # distanza tra nodi lungo x
Delta_y = H/Ny # distanza tra nodi lungo y
t_ele = 0.14 # spessore elementi (7+7 cm)
gamma_l = 25 # peso specifico (CLS = 25 kN/m^3)
N_ver = 400 # sforzo normale sul pannello
F_hor = 1340 # risultante orizzontale in cima al pannello
N_F_hor = int(np.ceil(0.3 / Delta_y)) + 1 # numero di punti in cui è suddivisa la forza orizzontale (distribuito su 30 cm)

nodes = []
for iy in range(Ny+1):
    for ix in range(Nx+1):
        nodes.append({'x': ix*Delta_x, 'y': iy*Delta_y, 'm': Delta_x*Delta_y*t_ele*gamma_l/9.81})

elements = []
for iy in range(Ny):
    for ix in range(Nx):
        i = (Nx+1)*(iy+0) + ix+0
        j = (Nx+1)*(iy+0) + ix+1
        k = (Nx+1)*(iy+1) + ix+1
        elements.append({'i':i, 'j':j, 'k':k, 'E':30e6, 'ni':0.2, 't':t_ele})
        i = (Nx+1)*(iy+0) + ix+0
        j = (Nx+1)*(iy+1) + ix+0
        k = (Nx+1)*(iy+1) + ix+1
        elements.append({'i':i, 'j':j, 'k':k, 'E':30e6, 'ni':0.2, 't':t_ele})

constraints = []
for i in range(Nx+1): constraints.append({'id':i, 'u':True, 'v':True}) # vincolati alla base

forces = []
for ix in range(Nx+1):
    id = (Ny)*(Nx+1)+ix
    forces.append({'id':id, 'f':[0,-N_ver/(Nx+1)]}) # pressione verticale sui nodi superiori
for iy in range(N_F_hor):
    forces.append({'id': (Nx+1)*(Ny-iy), 'f':[F_hor/N_F_hor,0]}) # forza orizzontale



for i,_ in enumerate(nodes): nodes[i]['id'] = i # aggiunge gli 'id' ai nodi (serve solo per facilitare agli umani la lettura degli elementi)

# importa il modello da un file JSON:
# with open("tri2D_model_2.json",'r') as f: data = json.load(f)
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



def add_element_stiffness(element):
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
        f[id*dof+0] = 0
    if constraint['v'] == True: # il grado di libertà verticale è vincolato
        for i in range(N*dof): K[id*dof+1,i] = 0; K[i,id*dof+1] = 0 # cancella tutte le righe e colonne associate al valore [id*dof+1,id*dof+1]
        K[id*dof+1,id*dof+1] = 1
        f[id*dof+1] = 0

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






# print('spostamenti:\n', u)
print(f"spostamento massimo: {np.max(u):.4f} m")
# ASSEMBLA DI NUOVO LA MATRICE K PER CALCOLARE LE FORZE....
K = np.zeros((N*dof, N*dof))
for e in elements: add_element_stiffness(e)
forze = K @ u
# print('forze:\n', forze)
for i in range(N):
    nodes[i]['fx'] = forze[i*dof].item()
    # if nodes[i]['x'] == 0: print(nodes[i]['y'], nodes[i]['fx'])



def compute_stress_criteria(elements):
    for e in elements:
        i, j, k = e['i'], e['j'], e['k']
        xi, yi = nodes[i]['x'], nodes[i]['y']
        xj, yj = nodes[j]['x'], nodes[j]['y']
        xk, yk = nodes[k]['x'], nodes[k]['y']
        
        E = e['E']
        ni = e['ni']
        t = e['t']
        
        A = abs((xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj)) / 2)
        b1, b2, b3 = yj - yk, yk - yi, yi - yj
        c1, c2, c3 = xk - xj, xi - xk, xj - xi
        B = (1/(2*A)) * np.array([
            [b1,  0, b2,  0, b3,  0],
            [ 0, c1,  0, c2,  0, c3],
            [c1, b1, c2, b2, c3, b3]
        ])
        D = (E/(1-ni*ni)) * np.array([[1, ni, 0], [ni, 1, 0], [0, 0, (1-ni)/2]])
        
        # vettore spostamenti locale
        ue = np.array([
            nodes[i]['u'], nodes[i]['v'],
            nodes[j]['u'], nodes[j]['v'],
            nodes[k]['u'], nodes[k]['v']
        ]).reshape(-1,1)
        
        sigma = D @ B @ ue
        sigma_x, sigma_y, tau_xy = sigma.flatten()
        
        # Von Mises
        sigma_vm = np.sqrt(sigma_x**2 - sigma_x*sigma_y + sigma_y**2 + 3*tau_xy**2)
        
        # Tresca (plane stress)
        sigma_1 = 0.5*(sigma_x + sigma_y) + np.sqrt(((sigma_x - sigma_y)/2)**2 + tau_xy**2) # sforzo principale 1
        sigma_2 = 0.5*(sigma_x + sigma_y) - np.sqrt(((sigma_x - sigma_y)/2)**2 + tau_xy**2) # sforzo principale 2
        tresca = max(abs(sigma_1 - sigma_2), abs(sigma_1), abs(sigma_2))
        
        e['sigma'] = {'sigma_x': sigma_x, 'sigma_y': sigma_y, 'tau_xy': tau_xy,
                      'von_mises': sigma_vm, 'tresca': tresca}


compute_stress_criteria(elements)

max_sigma = max(elements, key=lambda x: x['sigma']['von_mises'])['sigma']['von_mises']
min_sigma = min(elements, key=lambda x: x['sigma']['von_mises'])['sigma']['von_mises']

print(f"max(sigma_VM) = {max_sigma/1000:.2f} MPa")
# print(f"min(sigma_VM) = {min_sigma/1000:.2f} MPa")

# for idx, e in enumerate(elements):
    # print(f"Elemento {idx}: {e['sigma']}")







# exit()
# --------------------------------------------------------------------
# salva il modello:
model = {'model':'Tri 2D', 'units':units, 'author':'Andrea Marchi',  'nodes': nodes, 'elements': elements, 'constraints': constraints, 'forces': forces}
with open('lastra_model.json', 'w', encoding='utf-8') as f: json.dump(model, f, ensure_ascii=False, indent='\t')

# PLOT -------------------------------------------------------------------- per il debug??..........................

# indeformata:
# for e in elements: plt.plot([nodes[e['i']]['x'], nodes[e['j']]['x'], nodes[e['k']]['x'], nodes[e['i']]['x']], [nodes[e['i']]['y'], nodes[e['j']]['y'], nodes[e['k']]['y'], nodes[e['i']]['y']], '-o',color='lightgray')
# # STAMPARE ANCHE GLI ASSI LOCALI (almeno il vettore che va da "i" a "j")
# # deformata:
# if 'u' in nodes[0]: # controlla se ci sono gli spostamenti da mostrare
    # mul = 50 # scala che moltiplica gli spostamenti (multiplier)
    # for e in elements:
        # plt.plot([nodes[e['i']]['x']+mul*nodes[e['i']]['u'], nodes[e['j']]['x']+mul*nodes[e['j']]['u'], nodes[e['k']]['x']+mul*nodes[e['k']]['u'], nodes[e['i']]['x']+mul*nodes[e['i']]['u']], [nodes[e['i']]['y']+mul*nodes[e['i']]['v'], nodes[e['j']]['y']+mul*nodes[e['j']]['v'], nodes[e['k']]['y']+mul*nodes[e['k']]['v'], nodes[e['i']]['y']+mul*nodes[e['i']]['v']], '-ok')
    # for f in forces: # forze applicate
        # scale = 0.005
        # x = nodes[f['id']]['x']; u = nodes[f['id']]['u']
        # y = nodes[f['id']]['y']; v = nodes[f['id']]['v']
        # plt.arrow(x+mul*u, y+mul*v, f['f'][0]*scale, f['f'][1]*scale, head_width=2*scale, head_length=4*scale, color='red', length_includes_head=True)
        # plt.text(x + f['f'][0]*scale, y + f['f'][1]*scale, f"({f['f'][0]:.2f},{f['f'][1]:.2f})", color='red')
    # # for i,n in enumerate(nodes): # forze in tutti i nodi:
        # # scale = 0.2 # scala visualizzazione
        # # plt.arrow(n['x']+mul*n['u'], n['y']+mul*n['v'], forze[i*dof+0].item()*scale, forze[i*dof+1].item()*scale, head_width=0.02, head_length=0.04, color='red', length_includes_head=True)
        # # # plt.text(nodes[f['id']]['x'] + f['f'][0]*scale, nodes[f['id']]['y'] + f['f'][1]*scale, f"({f['f'][0]},{f['f'][1]})", color='red')
# plt.show()


# exit()

# # stampa la matrice K
# plt.imshow(K, cmap='viridis', origin='lower')
# for i in range(K.shape[0]):
    # for j in range(K.shape[1]):
        # plt.text(j, i, str(K[i, j]), ha='center', va='center', color='white') # mostra anche i valori
# plt.colorbar()  # aggiunge la scala dei valori
# plt.show()


# stampa gli sforzi di Von Mises

# --- estrai coordinate nodi ---
x = np.array([n['x'] for n in nodes])
y = np.array([n['y'] for n in nodes])

# --- lista triangoli (element connectivity) ---
tris = np.array([[e['i'], e['j'], e['k']] for e in elements])

# --- valori di von Mises per elemento ---
sigma_vm = np.array([e['sigma']['von_mises']/1000 for e in elements]) # in MPa

# triangolazione
triang = tri.Triangulation(x, y, tris)

# plot
plt.figure(figsize=(8, 6))
tpc = plt.tripcolor(triang,
                    sigma_vm,
                    shading='flat',    # 1 colore per elemento (CST)
                    cmap='gist_yarg')
# cmap = "inferno" (la più usata per FEA, viola → arancio → giallo → rosso)
# cmap = 'magma' (forse il migliore)
# cmap = "plasma"
# cmap = "viridis" (verde–blu–giallo)
# cmap = 'seismic'
# cmap = 'gist_rainbow' (arcobaleno)
# cmap = 'coolwarm' (da celeste a rosso)
# cmap = 'gist_yarg' (elegante, contrario di gist_gray)

# colorbar
plt.colorbar(tpc, label="Von Mises [MPa]")

plt.gca().set_aspect('equal')
plt.title("Mappa delle tensioni di Von Mises")
plt.xlabel("x [m]")
plt.ylabel("y [m]")

# plt.show()



# fattore globale (decidi quanto devono essere visibili le barre)
k = 0.000002   # scegli tu, 0.01–0.05 funziona bene in genere


for e in elements:
    i, j, k_node = e['i'], e['j'], e['k']

    # baricentro
    xc = (nodes[i]['x'] + nodes[j]['x'] + nodes[k_node]['x'])/3
    yc = (nodes[i]['y'] + nodes[j]['y'] + nodes[k_node]['y'])/3

    # sforzi dell'elemento
    sx = e['sigma']['sigma_x']
    sy = e['sigma']['sigma_y']
    txy = e['sigma']['tau_xy']

    # sforzi principali
    s1 = 0.5*(sx + sy) + np.sqrt(0.25*(sx - sy)**2 + txy**2)
    s2 = 0.5*(sx + sy) - np.sqrt(0.25*(sx - sy)**2 + txy**2)

    # direzione principale (angolo)
    theta = 0.5 * np.arctan2(2*txy, sx - sy)

    # versori delle direzioni principali
    dx1, dy1 = np.cos(theta), np.sin(theta)            # σ1
    dx2, dy2 = -np.sin(theta), np.cos(theta)           # σ2

    # lunghezze proporzionali al valore locale della tensione principale
    L1 = k * abs(s1)
    L2 = k * abs(s2)

    # barretta σ1 (nero)
    plt.plot([xc - L1*dx1, xc + L1*dx1],
             [yc - L1*dy1, yc + L1*dy1],
             color='black', linewidth=0.5)

    # barretta σ2 (bianco)
    plt.plot([xc - L2*dx2, xc + L2*dx2],
             [yc - L2*dy2, yc + L2*dy2],
             color='black', linewidth=0.5)

plt.gca().set_aspect('equal')
plt.title("Direzioni principali")
plt.xlabel("x")
plt.ylabel("y")
plt.show()






exit()
# ------------------------------------------------------------------------------------------------------------------------------------------------
# analisi modale:

# ricalcola la matrice delle rigidezze per l'analisi modale
K = np.zeros((N*dof, N*dof))
for e in elements: add_element_stiffness(e) # assembla la matrice di rigidezza K
for c in constraints: add_constraint(c) # aggiunge i vincoli (sempre alla matrice K)


# crea la matrice delle masse (concentrate nei nodi):
M = np.zeros((N*dof, N*dof)) # FARE CON LA MATRICE DELLE MASSE COMPATIBILE CON GLI ELEMENTI???
for i, node in enumerate(nodes):
    M[i*dof, i*dof] = node['m']      # grado di libertà orizzontale
    M[i*dof+1, i*dof+1] = node['m']  # grado di libertà verticale
# applica i vincoli alla matrice delle masse
for c in constraints:
    id = c['id']
    if c['u']: M[id*dof+0, id*dof+0] = 1 # oppure 0, a seconda di come vuoi trattare i vincoli
    if c['v']: M[id*dof+1, id*dof+1] = 1

# eliminiamo i gradi di libertà vincolati
free_dof = []
for i in range(N*dof):
    if K[i,i] != 1:  # i gradi di libertà vincolati in K li hai messi a 1
        free_dof.append(i)

K_reduced = K[np.ix_(free_dof, free_dof)]
M_reduced = M[np.ix_(free_dof, free_dof)]

# calcolo autovalori e autovettori
eigvals, eigvecs = eigh(K_reduced, M_reduced)

# autovalori positivi
omega = np.sqrt(eigvals)        # rad/s
frequencies = omega / (2*np.pi) # Hz
periods = 1 / frequencies       # s

for i in range(5):
    print(f"Modo {i+1}: T = {periods[i]:.5f} s")






# plotta la forma modale --------------------------------------------------------------------

mode = 0 # 0: primo modo
phi = eigvecs[:, mode]

# ricreiamo un array con tutti i gradi di libertà (vincolati = 0)
u_modal = np.zeros(N*dof)
u_modal[free_dof] = phi

# estrai spostamenti x e y
u_x = u_modal[0::2]
u_y = u_modal[1::2]

# scala per visualizzazione
scale = 0.2  # aumenta se il modo è troppo piccolo


plt.figure()
for e in elements: # modello indeformato
    x = [nodes[e['i']]['x'], nodes[e['j']]['x'], nodes[e['k']]['x'], nodes[e['i']]['x']]
    y = [nodes[e['i']]['y'], nodes[e['j']]['y'], nodes[e['k']]['y'], nodes[e['i']]['y']]
    plt.plot(x, y, '-o', color='lightgray')
    
for e in elements: # modo deformato
    x_def = [nodes[e['i']]['x'] + u_x[e['i']]*scale,
             nodes[e['j']]['x'] + u_x[e['j']]*scale,
             nodes[e['k']]['x'] + u_x[e['k']]*scale,
             nodes[e['i']]['x'] + u_x[e['i']]*scale]
    y_def = [nodes[e['i']]['y'] + u_y[e['i']]*scale,
             nodes[e['j']]['y'] + u_y[e['j']]*scale,
             nodes[e['k']]['y'] + u_y[e['k']]*scale,
             nodes[e['i']]['y'] + u_y[e['i']]*scale]
    plt.plot(x_def, y_def, '-ok')  # modo deformato

plt.title(f"Forma modale {mode+1}, T = {periods[mode]:.3f} s")
plt.axis('equal')
plt.show()
