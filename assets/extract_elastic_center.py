from abaqus import mdb
import math

model_name = mdb.models.keys()[3]
model = mdb.models[model_name]
assembly = model.rootAssembly
print('Running on model: %s' % model_name)

# Hardcode effective laminate Ex per member (psi)
E_skin = 8.3e6
E_spar = 20.0e6

t_skin = 0.02 + 0.11811 + 0.02
t_spar = 4 * 0.01

members = [
    ('Skin-1', t_skin, E_skin),
    ('FrontSpar-1', t_spar, E_spar),
    ('RearSpar-1',  t_spar, E_spar),
]

int_E_ds  = 0.0
int_Ey_ds = 0.0

for inst_name, t, E in members:
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
        int_E_ds  += E * t * ds
        int_Ey_ds += E * t * ym * ds

y_bar = int_Ey_ds / int_E_ds
print('Elastic center Y (chordwise): %.5f in' % y_bar)
print('As percent chord: %.1f%%' % (y_bar / 20.0 * 100.0))