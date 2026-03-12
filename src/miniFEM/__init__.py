# lo scopo principale di questo file è dire a Python "Questa cartella è un package Python".


# semantic versioning (semver.org)
VERSION = "0.2.0"

# Importa i moduli per renderli accessibili come sottomoduli (fem.tri2D, fem.truss2D)
from . import tri2D
from . import truss2D

# Esponi le funzioni principali di core.py direttamente nel namespace principale (fem.solve_system)
from .core import init_system, create_forces, add_constraint, solve_system



'''
ESEMPIO:

import miniFEM as fem

# 1. Dati in ingresso:
nodes = [...]
elements = [...]
constraints = [...]
forces = [...]

# 2. Inizializzazione:
K, f = fem.init_system(len(nodes))
fem.create_forces(f, forces)

# 3. Assemblaggio Matrice:
for e in elements:
    fem.tri2D.add_element_stiffness(K, e, nodes)  # Se è un triangolo CST
    # fem.truss2D.add_element_stiffness(K, e, nodes) # Se fosse una biella

# 4. Applicazione Vincoli:
for c in constraints:
    fem.add_constraint(K, f, c)

# 5. Risoluzione (con matrici sparse):
u = fem.solve_system(K, f)

print(u)
'''