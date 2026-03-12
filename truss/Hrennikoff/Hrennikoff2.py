import numpy as np
import matplotlib.pyplot as plt
import json

# aggiunge la cartella src al path (cartella dove sta miniFEM)
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../src'))

from miniFEM import truss2D as fem
import miniFEM as mfem

print('miniFEM version', mfem.VERSION)

# ----------------------------------------------------------------------------

# funzione che aggiunge gli elementi e i nodi (se serve) per approssimare un
# elemento finito quadrato per il continuo 2D
# il vettore "nodi" e' una lista di punti [(x1,y1),(x2,y2),...]
# aggiunge ai vettori nodes, elements
def aggiungi_elemento_Hrennikoff(nodi, t, E):
    # i nodi sono numerati in senso antiorario partendo da quello in basso a sinistra
    indici = [0,0,0,0] # indici dei quattro nodi dell'elemento
    # aggiunge i nodi se non sono presenti:
    for i,n in enumerate(nodi):
        aggiungi = True
        for ii,nn in enumerate(nodes):
            tmp = (nn['x'],nn['y'])
            if n == tmp:
                aggiungi = False
                indici[i] = ii
        if aggiungi:
            indici[i] = len(nodes)
            nodes.append({'x':n[0], 'y':n[1]})
    
    B = nodes[indici[1]]['x'] - nodes[indici[0]]['x'] # RICHIEDE CHE SIA NUMERATO PERFETTAMENTE
    H = nodes[indici[2]]['y'] - nodes[indici[1]]['y'] # E CHE SIA PERFETTAMENTE RETTANGOLARE
    D = np.sqrt(H**2 + B**2)
    k1 = 3/16*E*t/B*(3*B**2 - H**2)
    k2 = 3/16*E*t/H*(3*H**2 - B**2)
    k3 = 3/16*E*t*D**3/(B*H)
    elements.append({'i':indici[0], 'j':indici[1], 'E':E, 'A':k2/E}) # asta orizzontale
    elements.append({'i':indici[2], 'j':indici[3], 'E':E, 'A':k2/E}) # asta orizzontale
    elements.append({'i':indici[1], 'j':indici[2], 'E':E, 'A':k1/E}) # asta verticale
    elements.append({'i':indici[0], 'j':indici[3], 'E':E, 'A':k1/E}) # asta verticale
    elements.append({'i':indici[0], 'j':indici[2], 'E':E, 'A':k3/E}) # asta diagonale
    elements.append({'i':indici[1], 'j':indici[3], 'E':E, 'A':k3/E}) # asta diagonale


# MODELLO --------------------------------------------------------------------

units = "m, kN, kPa, s"


nodes = []
elements = []
Dx = 1
Dy = 1
for ix in range(0,5):
    for iy in range(0,5):
        aggiungi_elemento_Hrennikoff([
        ((ix+0)*Dx,(iy+0)*Dy),
        ((ix+1)*Dx,(iy+0)*Dy),
        ((ix+1)*Dx,(iy+1)*Dy),
        ((ix+0)*Dx,(iy+1)*Dy)], 
        0.1, 200e6)
# print(nodes)

# vincola il vertice in basso a sinistra:
constraints = [{'id':0, 'u':True, 'v':True}] # 'u':True vuol dire che lo spostamento lungo x è vincolato (come OpenSees)
# vincola tutti i nodi che sono a y=0 con carrelli che possono spostarsi lungo x:
for i,n in enumerate(nodes):
    if n['y'] == 0:
        constraints.append({'id':i, 'u':False, 'v':True})

forces = []
maxY = 0;
for n in nodes:
    if n['y'] > maxY: maxY = n['y']
for i,n in enumerate(nodes):
    if n['y'] == maxY:
        forces.append({'id':i, 'f':[0,-1]})

for i,_ in enumerate(nodes): nodes[i]['id'] = i # aggiunge gli 'id' ai nodi (serve solo per facilitare agli umani la lettura degli elementi)

# importa il modello da un file JSON:
# with open("truss2D_model.json",'r') as f: data = json.load(f)
# nodes = data['nodes']
# elements = data['elements']
# constraints = data['constraints']
# forces = data['forces']

# --------------------------------------------------------------------

dof = 2 # due spostamenti nel piano (no rotazioni)

# matrice di rigidezza del sistema
K = np.zeros((len(nodes)*dof, len(nodes)*dof)) # default: dtype=float64
# FARE CON LE MATRICI SPARSE (se il numero di nodi è eccessivo)


f = fem.create_forces(forces,nodes)
for e in elements: fem.add_element_stiffness(K, e, nodes) # assembla la matrice di rigidezza K
for c in constraints: fem.add_constraint(K, c, nodes) # aggiunge i vincoli (sempre alla matrice K)
u = np.linalg.solve(K,f) # trova gli spostamenti invertendo la matrice K (uguale a "u = np.linalg.inv(K) @ f", ma meno efficiente e meno numericamente stabile)


# salva gli spostamenti nei nodi:
for i,_ in enumerate(nodes):
    nodes[i]['u'] = u[i*dof+0].item()
    nodes[i]['v'] = u[i*dof+1].item()


# print('spostamenti:\n', u)
# ASSEMBLA DI NUOVO LA MATRICE K PER CALCOLARE LE FORZE....
K = np.zeros((len(nodes)*dof, len(nodes)*dof))
for e in elements: fem.add_element_stiffness(K, e, nodes)
forze = K @ u
# print('forze:\n', forze)

# dopo aver calcolato f = K @ u (servono u e v per i nodi)
element_forces = fem.compute_element_forces(elements, nodes)
# salva gli sforzi nelle aste:
for i,ef in enumerate(element_forces):
    elements[i]['N'] = ef['N']
    # print(f"Elemento {ef['element']['i']}-{ef['element']['j']}:  N = {ef['N']:.3f}   sigma = {ef['sigma']:.3f}")


# exit()
# --------------------------------------------------------------------
# salva il modello:
# model = {'model':'Truss 2D', 'units':units, 'author':'Andrea Marchi',  'nodes': nodes, 'elements': elements, 'constraints': constraints, 'forces': forces}
# with open('truss2D_model.json', 'w', encoding='utf-8') as f: json.dump(model, f, ensure_ascii=False, indent='\t')

# PLOT -------------------------------------------------------------------- per il debug??..........................

# indeformata:
for e in elements: plt.plot([nodes[e['i']]['x'], nodes[e['j']]['x']], [nodes[e['i']]['y'], nodes[e['j']]['y']], '-o', color='0.8')
# STAMPARE ANCHE GLI ASSI LOCALI (almeno il vettore che va da "i" a "j")
# deformata:
if 'u' in nodes[0]: # controlla se ci sono gli spostamenti da mostrare
    mul = 1000000 # scala che moltiplica gli spostamenti (multiplier)
    for e in elements:
        plt.plot([nodes[e['i']]['x']+mul*nodes[e['i']]['u'], nodes[e['j']]['x']+mul*nodes[e['j']]['u']],
                 [nodes[e['i']]['y']+mul*nodes[e['i']]['v'], nodes[e['j']]['y']+mul*nodes[e['j']]['v']], '-ok')
    for i,n in enumerate(nodes):
        scale = 0.8 # scala visualizzazione
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



# fem.plot_element_forces(nodes, element_forces)
# fem.plot_forces(nodes, forces)

# # VORREI FARE ANCHE CHE SALVA GLI SFORZI SULLE ASTE NEL JSON...
# # ANCHE LE FORZE NEI NODI


# plt.show()