#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 15:27:22 2026

@author: louislin
"""

from odbAccess import *
import math

# ── CONFIG ────────────────────────────────────────────────────────────────────
ODB_NAME  = 'Wing_Frequency.odb'
INST_NAME = 'SKIN-1'   # Abaqus uppercases; check with list(odb.rootAssembly.instances)
T         = 0.040      # skin thickness [in]  ← match SKIN_T

# ── EXTRACT ROOT CROSS-SECTION NODES (X ≈ 0) ─────────────────────────────────
odb  = openOdb(ODB_NAME)
inst = odb.rootAssembly.instances[INST_NAME]
root = [n for n in inst.nodes if n.coordinates[0] < 0.5]

# Sort by angle around profile centroid — handles closed airfoil correctly
yc = sum([n.coordinates[1] for n in root]) / len(root)
zc = sum([n.coordinates[2] for n in root]) / len(root)
root.sort(key=lambda n: math.atan2(n.coordinates[2] - zc,
                                    n.coordinates[1] - yc))

coords = [(n.coordinates[1], n.coordinates[2]) for n in root]
coords.append(coords[0])   # close loop

# ── SINGLE-PASS INTEGRALS ─────────────────────────────────────────────────────
C = s_y = s_z = s_y2 = s_z2 = shoelace = 0.0

for i in range(len(coords) - 1):
    y0, z0 = coords[i];  y1, z1 = coords[i + 1]
    ds = math.sqrt((y1 - y0)**2 + (z1 - z0)**2)
    ym = 0.5 * (y0 + y1);  zm = 0.5 * (z0 + z1)
    C        += ds
    s_y      += ym * ds
    s_z      += zm * ds
    s_y2     += ym**2 * ds
    s_z2     += zm**2 * ds
    shoelace += y0 * z1 - y1 * z0   # shoelace for enclosed area

y_bar = s_y / C
z_bar = s_z / C
A_enc = 0.5 * abs(shoelace)

# ── SECTION PROPERTIES ────────────────────────────────────────────────────────
A   = T * C
Ixx = T * (s_z2 - C * z_bar**2)   # about centroidal chord axis
Iyy = T * (s_y2 - C * y_bar**2)   # about centroidal thickness axis
Ip  = Ixx + Iyy
J   = 4.0 * A_enc**2 * T / C      # Bredt-Batho (uniform t)

print('Nodes  : %d' % len(root))
print('Perim  : %.4f in'   % C)
print('A_enc  : %.4f in^2' % A_enc)
print('A      : %.5f in^2' % A)
print('Ixx    : %.5f in^4' % Ixx)
print('Iyy    : %.5f in^4' % Iyy)
print('J      : %.5f in^4' % J)
print('Ip     : %.5f in^4' % Ip)

odb.close()