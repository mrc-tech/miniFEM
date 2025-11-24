import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.linalg import eigh

# -------------------
# Parametri serbatoio
# -------------------
R = 43       # raggio serbatoio [m]
H = 20       # altezza liquido [m]
Nr = 20        # numero di cerchi concentrici
Ntheta = 64   # numero di punti per cerchio
rho = 857.0  # densità [kg/m^3]
g = 9.81      # gravità [m/s^2]
t_ele = 1   # spessore unitario

# -------------------
# Generazione nodi in coordinate polari e conversione cartesiana
# -------------------
nodes = []
for ir in range(Nr+1):
    r = R * ir / Nr
    ntheta = max(1, int(Ntheta * ir / Nr))  # centro ha 1 nodo
    for itheta in range(ntheta):
        theta = 2*np.pi * itheta / ntheta
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        nodes.append({'x': x, 'y': y, 'r': r, 'theta': theta})

N = len(nodes)

# -------------------
# Triangolazione con matplotlib.tri
# -------------------
x = np.array([n['x'] for n in nodes])
y = np.array([n['y'] for n in nodes])
triang = tri.Triangulation(x, y)
elements = [{'i': t[0], 'j': t[1], 'k': t[2]} for t in triang.triangles]

# -------------------
# Inizializzazione matrici M e K
# -------------------
M = np.zeros((N, N))
K = np.zeros((N, N))

for e in elements:
    i, j, k = e['i'], e['j'], e['k']
    xi, yi = x[i], y[i]
    xj, yj = x[j], y[j]
    xk, yk = x[k], y[k]
    
    # area triangolo
    A = abs((xi*(yj-yk) + xj*(yk-yi) + xk*(yi-yj)) / 2)
    
    # matrice massa lumped
    m_ele = rho * t_ele * A / 3
    for n in [i,j,k]:
        M[n,n] += m_ele
    
    # matrice rigidezza CST lineare (effetto gravità)
    b = np.array([yj-yk, yk-yi, yi-yj])
    c = np.array([xk-xj, xi-xk, xj-xi])
    ke = np.zeros((3,3))
    for a in range(3):
        for b_ in range(3):
            ke[a,b_] = (b[a]*b[b_] + c[a]*c[b_]) / (4*A) * rho * g * t_ele
    dof_map = [i,j,k]
    for a in range(3):
        for b_ in range(3):
            K[dof_map[a], dof_map[b_]] += ke[a,b_]

# -------------------
# Vincoli
# -------------------
fixed_nodes = []
for idx, n in enumerate(nodes):
    if np.isclose(n['r'], R, atol=1e-6):  # parete verticale
        fixed_nodes.append(idx)

free_dof = [i for i in range(N) if i not in fixed_nodes]

K_reduced = K[np.ix_(free_dof, free_dof)]
M_reduced = M[np.ix_(free_dof, free_dof)]

# -------------------
# Analisi modale
# -------------------
eigvals, eigvecs = eigh(K_reduced, M_reduced)
omega = np.sqrt(eigvals)
frequencies = omega/(2*np.pi)
periods = 1/frequencies

print("Prime 5 frequenze e periodi (slosh in serbatoio circolare):")
for i in range(5):
    print(f"Modo {i+1}: f = {frequencies[i]:.3f} Hz, T = {periods[i]:.3f} s")

# -------------------
# Plot prima forma modale
# -------------------
mode = 0
phi = eigvecs[:, mode]
u_modal = np.zeros(N)
u_modal[free_dof] = phi

scale = 10
plt.figure()
for tri_idx in triang.triangles:
    x_tri = x[tri_idx]
    y_tri = y[tri_idx]
    plt.plot(np.append(x_tri, x_tri[0]), np.append(y_tri, y_tri[0]), '-o', color='lightgray')
for tri_idx in triang.triangles:  
    x_tri = x[tri_idx]
    y_tri = y[tri_idx]
    y_def = y_tri + u_modal[tri_idx]*scale
    plt.plot(np.append(x_tri, x_tri[0]), np.append(y_def, y_def[0]), '-ok')

plt.title(f"Prima forma modale sloshing (serbatoio circolare), T = {periods[mode]:.3f} s")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.axis('equal')
plt.show()
