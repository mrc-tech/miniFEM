import numpy as np
from scipy.sparse import lil_matrix  # matrici sparse LIL sono piu' semplici per aggiungere elementi
from scipy.sparse.linalg import spsolve

dof = 2 # Gradi di libertà per nodo

def init_system(num_nodes):
    ''' Inizializza la matrice di rigidezza sparsa K e il vettore delle forze f. '''
    size = num_nodes * dof
    # lil_matrix è la più veloce per inserire/modificare elementi uno ad uno
    K = lil_matrix((size, size), dtype=np.float64)
    f = np.zeros((size, 1), dtype=np.float64)
    return K, f

def create_forces(f, forces):
    ''' Somma le forze nodali al vettore globale f. '''
    for x in forces:
        f[x['id']*dof+0, 0] += x['f'][0]
        f[x['id']*dof+1, 0] += x['f'][1]

def add_constraint(K, f, constraint):
    '''
        Applica i vincoli:
        1. Azzera righe e colonne in K, mettendo 1 sulla diagonale.
        2. Azzera il corrispondente valore nel vettore delle forze f.
    '''
    
    if constraint.get('u', False): # Vincolo orizzontale
        idx = constraint['id'] * dof + 0
        K[idx, :] = 0
        K[:, idx] = 0 # <--- CRITICO !!!!
        '''
            Hai scelto (giustamente) la lil_matrix di SciPy. LIL sta per List of Lists, ed è una struttura dati salvata per righe.
            Quando scrivi K[idx, :] = 0, Python cancella una riga in una frazione di millisecondo.
            Ma quando scrivi K[:, idx] = 0 (cancellare una colonna), costringi SciPy a fare un ciclo su tutte le righe della matrice per cercare ed eliminare quel singolo elemento. Su una matrice 10.000 x 10.000, azzerare le colonne per applicare i vincoli farà letteralmente "congelare" il tuo computer per diversi secondi o minuti.
        '''
        K[idx, idx] = 1.0
        f[idx, 0] = 0.0 # Il valore del vettore f è ora nullo!
        
    if constraint.get('v', False): # Vincolo verticale
        idx = constraint['id'] * dof + 1
        K[idx, :] = 0
        K[:, idx] = 0
        K[idx, idx] = 1.0
        f[idx, 0] = 0.0

def add_constraint_penalty(K, f, constraint):
    '''
        Applica i vincoli usando il Metodo della Penalità per alte prestazioni
        con matrici sparse.
    '''
    penalty = 1e15  # Numero enorme
    
    if constraint.get('u', False): # Vincolo orizzontale
        idx = constraint['id'] * dof + 0
        K[idx, idx] += penalty
        # Se volessimo imporre uno spostamento (cedimento) di 5m, scriveremmo: 
        # f[idx, 0] += 5.0 * penalty. 
        # Siccome il vincolo è fisso (0 spostamento), lasciamo f inalterato o lo poniamo a 0
        f[idx, 0] = 0.0 
        
    if constraint.get('v', False): # Vincolo verticale
        idx = constraint['id'] * dof + 1
        K[idx, idx] += penalty
        f[idx, 0] = 0.0

def solve_system(K, f):
    ''' Converte K nel formato CSR (ottimizzato per i calcoli) e risolve il sistema. '''
    K_csr = K.tocsr()
    u = spsolve(K_csr, f)
    return u.reshape(-1, 1) # Restituisce un vettore colonna