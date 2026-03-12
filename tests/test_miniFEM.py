# aggiunge la cartella src al path (cartella dove sta miniFEM)
import sys; import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))

import miniFEM as fem

# 1. Definizione modello:
nodes = [{'x':0, 'y':0}, {'x':1, 'y':0}, {'x':1, 'y':1}]
elements = [{'i':0, 'j':1, 'k':2, 'E':210e3, 'ni':0.3, 't':0.1}] # (tri2D)
# elements = [{'i':0, 'j':2, 'E':210e3, 'A':0.1}, {'i':1, 'j':2, 'E':210e3, 'A':0.1}] # (truss2D)
constraints = [{'id':0, 'u':True, 'v':True}, {'id':1, 'u':True, 'v':True}]
forces = [{'id':2, 'f':[1,0]}]

# 2. Inizializzazione :
K, f = fem.init_system(len(nodes))
fem.create_forces(f, forces)

# 3. Assemblaggio Matrice:
for e in elements:
    fem.tri2D.add_element_stiffness(K, e, nodes) # Se è un triangolo CST
    # fem.truss2D.add_element_stiffness(K, e, nodes) # Se fosse una biella

# 4. Applicazione Vincoli:
for c in constraints:
    fem.add_constraint(K, f, c)

# 5. Risoluzione (con matrici sparse):
u = fem.solve_system(K, f)

# 6. Salvataggio risultati nei nodi:
for i, _ in enumerate(nodes):
    nodes[i]['u'] = u[i*2 + 0].item()
    nodes[i]['v'] = u[i*2 + 1].item()

# 7. output:
print(u)
fem.tri2D.plot_displacement(nodes, elements, scale=1000) # tri2D
# fem.truss2D.plot_element_forces(nodes, fem.truss2D.compute_element_forces(elements, nodes)) # truss2D
