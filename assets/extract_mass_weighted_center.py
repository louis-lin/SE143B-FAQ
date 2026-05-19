from abaqus import mdb
import math

model_name = mdb.models.keys()[3] # CHANGE IF YOU HAVE MORE THAN 1 MODEL, 0 IS FIRST ONE.
model = mdb.models[model_name]
assembly = model.rootAssembly
print('Running on model: %s' % model_name)

# INPUTS 
# Hardcode effective laminate Ex per member (psi)
rho_skin = 0.056   # CF woven
rho_spar = 0.056   # CF UD

t_skin = 0.02 + 0.11811 + 0.02 # THICKNESS OF SKIN INSTANCE
t_spar = 4 * 0.01 ## THICKNESS OF SPAR INSTANCE

members = [ # INSTANCE NAME, THICKNESS OF SKIN, MODULUS OF SKIN
    ('Skin-1', t_skin, E_skin), 
    ('FrontSpar-1', t_spar, E_spar),
    ('RearSpar-1',  t_spar, E_spar),
]

# CALCULATIONS ....  D
int_rho_ds  = 0.0
int_rho_y_ds = 0.0

for inst_name, t, rho in members:
    inst = assembly.instances[inst_name]
    root_nodes = [n for n in inst.nodes if n.coordinates[0] < 0.5]
    root_nodes = sorted(root_nodes, key=lambda n: n.coordinates[1])
    for i in range(len(root_nodes) - 1):
        y0 = root_nodes[i].coordinates[1]
        y1 = root_nodes[i+1].coordinates[1]
        z0 = root_nodes[i].coordinates[2]
        z1 = root_nodes[i+1].coordinates[2]
        ds = math.sqrt((y1-y0)**2 + (z1-z0)**2)
        ym = 0.5*(y0 + y1)
        int_rho_ds   += rho * t * ds
        int_rho_y_ds += rho * t * ym * ds

y_cg = int_rho_y_ds / int_rho_ds

print('Section CG Y (chordwise): %.5f in' % y_cg)
print('As percent chord: %.1f%%' % (y_cg / 20.0 * 100.0))