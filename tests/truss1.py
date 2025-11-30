import numpy as np
import matplotlib.pyplot as plt
import json

# aggiunge la cartella src al path (cartella dove sta miniFEM)
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))

from miniFEM import truss2D as fem
import miniFEM as mfem

print(mfem.VERSION)

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


print('spostamenti:\n', u)
# ASSEMBLA DI NUOVO LA MATRICE K PER CALCOLARE LE FORZE....
K = np.zeros((len(nodes)*dof, len(nodes)*dof))
for e in elements: fem.add_element_stiffness(K, e, nodes)
forze = K @ u
print('forze:\n', forze)

# dopo aver calcolato f = K @ u (servono u e v per i nodi)
element_forces = fem.compute_element_forces(elements, nodes)
# salva gli sforzi nelle aste:
for i,ef in enumerate(element_forces):
    elements[i]['N'] = ef['N']
    print(f"Elemento {ef['element']['i']}-{ef['element']['j']}:  N = {ef['N']:.3f}   sigma = {ef['sigma']:.3f}")


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



fem.plot_element_forces(nodes, element_forces)
fem.plot_forces(nodes, forces)

# VORREI FARE ANCHE CHE SALVA GLI SFORZI SULLE ASTE NEL JSON...
# ANCHE LE FORZE NEI NODI


plt.show()