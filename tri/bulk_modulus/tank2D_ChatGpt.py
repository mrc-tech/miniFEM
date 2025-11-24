import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

# -------------------
# Parametri serbatoio
# -------------------
H = 20       # altezza
B = 86       # base
Nx = 86       # nodi in x
Ny = 20       # nodi in y
rho = 857.0  # kg/m^3
g = 9.81
t_ele = 1.0   # spessore unitario

Delta_x = B/Nx
Delta_y = H/Ny

# -------------------
# Nodi
# -------------------
nodes = []
for iy in range(Ny+1):
    for ix in range(Nx+1):
        nodes.append({'x': ix*Delta_x, 'y': iy*Delta_y})

N = len(nodes)
dof = 1  # pressione/altezza

# -------------------
# Elementi triangolari
# -------------------
elements = []
for iy in range(Ny):
    for ix in range(Nx):
        i = (Nx+1)*(iy+0) + ix+0
        j = (Nx+1)*(iy+0) + ix+1
        k = (Nx+1)*(iy+1) + ix+1
        elements.append({'i':i, 'j':j, 'k':k})
        i = (Nx+1)*(iy+0) + ix+0
        j = (Nx+1)*(iy+1) + ix+0
        k = (Nx+1)*(iy+1) + ix+1
        elements.append({'i':i, 'j':j, 'k':k})

# -------------------
# Matrici M e K
# -------------------
M = np.zeros((N*dof, N*dof))
K = np.zeros((N*dof, N*dof))

for e in elements:
    i, j, k = e['i'], e['j'], e['k']
    xi, yi = nodes[i]['x'], nodes[i]['y']
    xj, yj = nodes[j]['x'], nodes[j]['y']
    xk, yk = nodes[k]['x'], nodes[k]['y']
    
    # area triangolo
    A = abs((xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj)) / 2)
    
    # ---- matrice di massa concentrata ai nodi ----
    m_ele = rho * t_ele * A / 3
    for n in [i,j,k]:
        M[n,n] += m_ele
    
    # ---- matrice di rigidezza CST lineare ----
    # coefficienti b,c
    b = np.array([yj-yk, yk-yi, yi-yj])
    c = np.array([xk-xj, xi-xk, xj-xi])
    ke = np.zeros((3,3))
    for a in range(3):
        for b_ in range(3):
            ke[a,b_] = (b[a]*b[b_] + c[a]*c[b_]) / (4*A) * rho * g * t_ele
    # assemblaggio
    dof_map = [i,j,k]
    for a in range(3):
        for b_ in range(3):
            K[dof_map[a], dof_map[b_]] += ke[a,b_]

# -------------------
# Vincoli
# -------------------
fixed_nodes = []
for i, node in enumerate(nodes):
    if node['y'] == 0.0:  # fondo
        fixed_nodes.append(i)
    # if node['x'] == 0.0 or node['x'] == B:  # pareti verticali
        # fixed_nodes.append(i)

free_dof = [i for i in range(N*dof) if i not in fixed_nodes]

K_reduced = K[np.ix_(free_dof, free_dof)]
M_reduced = M[np.ix_(free_dof, free_dof)]

# -------------------
# Analisi modale
# -------------------
eigvals, eigvecs = eigh(K_reduced, M_reduced)
omega = np.sqrt(eigvals)
frequencies = omega/(2*np.pi)
periods = 1/frequencies

print("Prime 5 frequenze e periodi (slosh):")
for i in range(5):
    print(f"Modo {i+1}: f = {frequencies[i]:.3f} Hz, T = {periods[i]:.3f} s")

# -------------------
# Plot prima forma modale
# -------------------
mode = 1
phi = eigvecs[:, mode]
u_modal = np.zeros(N*dof)
u_modal[free_dof] = phi

scale = 3000
plt.figure()
for e in elements:
    x = [nodes[e['i']]['x'], nodes[e['j']]['x'], nodes[e['k']]['x'], nodes[e['i']]['x']]
    y = [nodes[e['i']]['y'], nodes[e['j']]['y'], nodes[e['k']]['y'], nodes[e['i']]['y']]
    plt.plot(x, y, '-o', color='lightgray')
for e in elements:
    y_def = [nodes[e['i']]['y'] + u_modal[e['i']]*scale,
             nodes[e['j']]['y'] + u_modal[e['j']]*scale,
             nodes[e['k']]['y'] + u_modal[e['k']]*scale,
             nodes[e['i']]['y'] + u_modal[e['i']]*scale]
    plt.plot([nodes[e['i']]['x'], nodes[e['j']]['x'], nodes[e['k']]['x'], nodes[e['i']]['x']], y_def, '-ok')

plt.title(f"Prima forma modale sloshing, T = {periods[mode]:.3f} s")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.axis('equal')
plt.show()
