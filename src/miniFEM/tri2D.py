import numpy as np
import json  # per importare ed esportare i dati del modello
import matplotlib.pyplot as plt  # per plottare il modello


# Le variabili globali vengono viste solo dalle funzioni del modulo
dof = 2 # gradi di liberta' per ogni nodo


# aggiunge la matrice di rigidezza dell'elemento alla matrice di rigidezza K (globale)
# def add_element_stiffness(K: np.ndarray, element: dict[str,float], nodes: dict[str,float]) -> None: # funziona per Python >= 3.9
def add_element_stiffness(K: np.ndarray, element, nodes) -> None:
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



def plot_displacement(nodes, elements, scale=100):
    # indeformata:
    for e in elements: plt.plot([nodes[e['i']]['x'], nodes[e['j']]['x'], nodes[e['k']]['x'], nodes[e['i']]['x']], [nodes[e['i']]['y'], nodes[e['j']]['y'], nodes[e['k']]['y'], nodes[e['i']]['y']], '-ok')
    # STAMPARE ANCHE GLI ASSI LOCALI (almeno il vettore che va da "i" a "j")
    # deformata:
    if 'u' in nodes[0]: # controlla se ci sono gli spostamenti da mostrare
        mul = scale # scala che moltiplica gli spostamenti (multiplier)
        for e in elements:
            plt.plot([nodes[e['i']]['x']+mul*nodes[e['i']]['u'], nodes[e['j']]['x']+mul*nodes[e['j']]['u'], nodes[e['k']]['x']+mul*nodes[e['k']]['u'], nodes[e['i']]['x']+mul*nodes[e['i']]['u']], [nodes[e['i']]['y']+mul*nodes[e['i']]['v'], nodes[e['j']]['y']+mul*nodes[e['j']]['v'], nodes[e['k']]['y']+mul*nodes[e['k']]['v'], nodes[e['i']]['y']+mul*nodes[e['i']]['v']], '-or')
        # for i,n in enumerate(nodes):
            # scala = 0.2 # scala visualizzazione
            # plt.arrow(n['x']+mul*n['u'], n['y']+mul*n['v'], forze[i*dof+0].item()*scale, forze[i*dof+1].item()*scale, head_width=0.02, head_length=0.04, color='red', length_includes_head=True)
            # # plt.text(nodes[f['id']]['x'] + f['f'][0]*scale, nodes[f['id']]['y'] + f['f'][1]*scale, f"({f['f'][0]},{f['f'][1]})", color='red')
    plt.show()

















# units = "m, kN, kPa, Mg, s"

# nodes = [
    # {'x':0, 'y':0},
    # {'x':1, 'y':0},
    # {'x':1, 'y':1}]

# elements = [
    # {'i':0, 'j':1, 'k':2, 'E':210e3, 'ni':0.3, 't':0.1}]

# constraints = [
    # {'id':0, 'u':True, 'v':True}, # 'u':True vuol dire che lo spostamento lungo x è vincolato (come OpenSees)
    # {'id':1, 'u':True, 'v':True}]

# forces = [
    # {'id':2, 'f':[1,0]}]

# for i,_ in enumerate(nodes): nodes[i]['id'] = i # aggiunge gli 'id' ai nodi (serve solo per facilitare agli umani la lettura degli elementi)

# # importa il modello da un file JSON:
# # with open("tri2D_model.json",'r') as f: data = json.load(f)
# # nodes = data['nodes']
# # elements = data['elements']
# # constraints = data['constraints']
# # forces = data['forces']

# # --------------------------------------------------------------------

# N = len(nodes) # numero di nodi
# dof = 2 # due spostamenti nel piano (no rotazioni)


# # print('nodes:'); print(nodes)
# # print('elements:'); print(elements)

# # matrice di rigidezza del sistema
# K = np.zeros((N*dof, N*dof)) # default: dtype=float64
# # FARE CON LE MATRICI SPARSE (se il numero di nodi è eccessivo)



# def add_element_stiffness(element):
    # """
    # Aggiunge la matrice di rigidezza dell'elemento triangolare 2D (CST) alla matrice globale K
    # """
    # i, j, k = element['i'], element['j'], element['k'] # nodi dell'elemento
    # xi, yi = nodes[i]['x'], nodes[i]['y']
    # xj, yj = nodes[j]['x'], nodes[j]['y']
    # xk, yk = nodes[k]['x'], nodes[k]['y']
    
    # E = element['E']
    # ni = element['ni']
    # t = element['t'] # spessore dell'elemento
    
    # A = abs((xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj)) / 2) # area del triangolo (sempre positiva)
    # if A == 0: raise ValueError("Triangolo degenerato con area zero!")
    
    # # --- coefficienti b e c ---
    # b1, b2, b3 = yj - yk, yk - yi, yi - yj
    # c1, c2, c3 = xk - xj, xi - xk, xj - xi
    # B = (1/(2*A)) * np.array([
        # [b1,  0, b2,  0, b3,  0],
        # [ 0, c1,  0, c2,  0, c3],
        # [c1, b1, c2, b2, c3, b3]
    # ])
    
    # D = (E/(1-ni*ni)) * np.array([[1, ni, 0], [ni, 1, 0], [0, 0, (1-ni)/2]]) # plane stress
    
    # ke = B.T @ D @ B * t * A # matrice di rigidezza dell'elemento (ke = B^T D B * A*t) ("A @ B" è un modo semplice per "np.matmul(A,B)")
    
    # # assemblaggio nella matrice globale K
    # dof_map = [i*dof, i*dof+1, j*dof, j*dof+1, k*dof, k*dof+1]
    # for r in range(6):
        # for c in range(6):
            # K[dof_map[r], dof_map[c]] += ke[r, c]

# def add_constraint(constraint): # aggiunge il vincolo
    # id = constraint['id']
    # if constraint['u'] == True: # il grado di libertà orizzontale è vincolato
        # for i in range(N*dof): K[id*dof+0,i] = 0; K[i,id*dof+0] = 0 # cancella tutte le righe e colonne associate al valore [id*dof+0,id*dof+0]
        # K[id*dof+0,id*dof+0] = 1
    # if constraint['v'] == True: # il grado di libertà verticale è vincolato
        # for i in range(N*dof): K[id*dof+1,i] = 0; K[i,id*dof+1] = 0 # cancella tutte le righe e colonne associate al valore [id*dof+1,id*dof+1]
        # K[id*dof+1,id*dof+1] = 1
    # # DOVREBBE CONTROLLARE CHE ANCHE I RELATIVI VALORI DEL VETTORE f SIANO NULLI!!!!!

# def create_forces(forces): # crea il vettore delle forze applicate
    # f = np.zeros((N*dof,1))
    # for x in forces: f[x['id']*dof+0], f[x['id']*dof+1] = x['f']
    # return f


# f = create_forces(forces)
# for e in elements: add_element_stiffness(e) # assembla la matrice di rigidezza K
# for c in constraints: add_constraint(c) # aggiunge i vincoli (sempre alla matrice K)
# u = np.linalg.solve(K,f) # trova gli spostamenti invertendo la matrice K (uguale a "u = np.linalg.inv(K) @ f", ma meno efficiente e meno numericamente stabile)


# # salva gli spostamenti nei nodi:
# for i,_ in enumerate(nodes):
    # nodes[i]['u'] = u[i*dof+0].item()
    # nodes[i]['v'] = u[i*dof+1].item()


# print('spostamenti:\n', u)
# # ASSEMBLA DI NUOVO LA MATRICE K PER CALCOLARE LE FORZE....
# K = np.zeros((N*dof, N*dof))
# for e in elements: add_element_stiffness(e)
# forze = K @ u
# print('forze:\n', forze)

# # CALCOLARE LE TENSIONI NEGLI ELEMENTI........



# def compute_stress_criteria(elements):
    # for e in elements:
        # i, j, k = e['i'], e['j'], e['k']
        # xi, yi = nodes[i]['x'], nodes[i]['y']
        # xj, yj = nodes[j]['x'], nodes[j]['y']
        # xk, yk = nodes[k]['x'], nodes[k]['y']
        
        # E = e['E']
        # ni = e['ni']
        # t = e['t']
        
        # A = abs((xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj)) / 2)
        # b1, b2, b3 = yj - yk, yk - yi, yi - yj
        # c1, c2, c3 = xk - xj, xi - xk, xj - xi
        # B = (1/(2*A)) * np.array([
            # [b1,  0, b2,  0, b3,  0],
            # [ 0, c1,  0, c2,  0, c3],
            # [c1, b1, c2, b2, c3, b3]
        # ])
        # D = (E/(1-ni*ni)) * np.array([[1, ni, 0], [ni, 1, 0], [0, 0, (1-ni)/2]])
        
        # # vettore spostamenti locale
        # ue = np.array([
            # nodes[i]['u'], nodes[i]['v'],
            # nodes[j]['u'], nodes[j]['v'],
            # nodes[k]['u'], nodes[k]['v']
        # ]).reshape(-1,1)
        
        # sigma = D @ B @ ue
        # sigma_x, sigma_y, tau_xy = sigma.flatten()
        
        # # Von Mises
        # sigma_vm = np.sqrt(sigma_x**2 - sigma_x*sigma_y + sigma_y**2 + 3*tau_xy**2)
        
        # # Tresca (plane stress)
        # sigma_1 = 0.5*(sigma_x + sigma_y) + np.sqrt(((sigma_x - sigma_y)/2)**2 + tau_xy**2) # sforzo principale 1
        # sigma_2 = 0.5*(sigma_x + sigma_y) - np.sqrt(((sigma_x - sigma_y)/2)**2 + tau_xy**2) # sforzo principale 2
        # tresca = max(abs(sigma_1 - sigma_2), abs(sigma_1), abs(sigma_2))
        
        # e['sigma'] = {'sigma_x': sigma_x, 'sigma_y': sigma_y, 'tau_xy': tau_xy,
                      # 'von_mises': sigma_vm, 'tresca': tresca}


# compute_stress_criteria(elements)

# for idx, e in enumerate(elements):
    # print(f"Elemento {idx}: {e['sigma']}")







# # exit()
# # --------------------------------------------------------------------
# # salva il modello:
# model = {'model':'Tri 2D', 'units':units, 'author':'Andrea Marchi',  'nodes': nodes, 'elements': elements, 'constraints': constraints, 'forces': forces}
# with open('tri2D_model.json', 'w', encoding='utf-8') as f: json.dump(model, f, ensure_ascii=False, indent='\t')

# # PLOT -------------------------------------------------------------------- per il debug??..........................



