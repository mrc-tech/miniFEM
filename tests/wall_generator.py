import json # per importare ed esportare i dati del modello
import matplotlib.pyplot as plt # PER IL DEBUG???


units = "m, kN, kPa, Mg, s"

B = 1 # base
H = 3 # altezza
Nx = 10 # numero di suddivisioni lungo x
Ny = 30 # numero di suddivisioni lungo y
Delta_x = B/Nx # distanza tra nodi lungo x
Delta_y = H/Ny # distanza tra nodi lungo y
E = 210e6 # Modulo elasticita'
A = 0.1 # Area aste

nodes = []
for iy in range(Ny+1):
    for ix in range(Nx+1):
        nodes.append({'x': ix*Delta_x, 'y': iy*Delta_y})

elements = []
for iy in range(Ny):
    for ix in range(Nx):
        i = (Nx+1)*(iy+0) + ix+0
        j = (Nx+1)*(iy+0) + ix+1
        elements.append({'i':i, 'j':j, 'E':E, 'A':A}) # orizzontale
        i = (Nx+1)*(iy+0) + ix+0
        j = (Nx+1)*(iy+1) + ix+0
        elements.append({'i':i, 'j':j, 'E':E, 'A':A}) # verticale
        i = (Nx+1)*(iy+0) + ix+1
        j = (Nx+1)*(iy+1) + ix+1
        elements.append({'i':i, 'j':j, 'E':E, 'A':A}) # verticale
        i = (Nx+1)*(iy+0) + ix+0
        j = (Nx+1)*(iy+1) + ix+1
        elements.append({'i':i, 'j':j, 'E':E, 'A':A}) # diagonale
        i = (Nx+1)*(iy+0) + ix+1
        j = (Nx+1)*(iy+1) + ix+0
        elements.append({'i':i, 'j':j, 'E':E, 'A':A}) # diagonale
for ix in range(Nx): # elementi superiori (orizzontali)
    i = (Nx+1)*(Ny) + ix+0
    j = (Nx+1)*(Ny) + ix+1
    elements.append({'i':i, 'j':j, 'E':E, 'A':A})

constraints = []
for i in range(Nx+1): constraints.append({'id':i, 'u':True, 'v':True}) # incastri alla base



# forces = []
forces = [
    {'id': (Nx+1)*(Ny), 'f':[1,0]}]


for i,_ in enumerate(nodes): nodes[i]['id'] = i # aggiunge gli 'id' ai nodi (serve solo per facilitare agli umani la lettura degli elementi)


# plot
fig, ax = plt.subplots()
for e in elements:
    ax.plot([nodes[e['i']]['x'], nodes[e['j']]['x']], [nodes[e['i']]['y'], nodes[e['j']]['y']], '-ok')
ax.set_aspect('equal')
plt.show()





# salva il modello:
model = {'model':'Truss 2D', 'units':units, 'author':'Andrea Marchi',  'nodes': nodes, 'elements': elements, 'constraints': constraints, 'forces': forces}
with open('wall.json', 'w', encoding='utf-8') as f: json.dump(model, f, ensure_ascii=False, indent='\t')
