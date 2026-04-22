# -*- coding: utf-8 -*-
"""
wing_abaqus.py
X = span (0=root, SPAN=tip), Y = chord (0=LE, CHORD=TE), Z = flap/thickness
Single closed skin shell, box-spar structure, N evenly-spaced ribs.
"""

from abaqus import *
from abaqusConstants import *
import regionToolset, mesh

# ── PARAMETERS ────────────────────────────────────────────────────────────────
CHORD    = 20.0
SPAN     = 84.0
YSPAR1   = 5.0    # front spar: 25% chord
YSPAR2   = 14.0   # rear spar:  70% chord
N_RIBS   = 8      # total ribs including root (X=0) and tip (X=SPAN)

E        = 10.0e6
NU       = 0.33
RHO      = 0.1 / 386.4
SKIN_T   = 0.040
SPAR_T   = 0.060
RIB_T    = 0.040
N_EIGEN  = 10
MESH_SIZE = 1.0

MODEL_NAME = 'Wing'
JOB_NAME   = 'Wing_Frequency'

RIB_X = [i * SPAN / (N_RIBS - 1) for i in range(N_RIBS)]

# ── AIRFOIL (Y=chord, Z=thickness) ───────────────────────────────────────────
_UPPER = [
    (0.00000, 0.00000),(0.00066, 0.00581),(0.00296, 0.01331),(0.00688, 0.02129),
    (0.01236, 0.02958),(0.01937, 0.03807),(0.02789, 0.04665),(0.03787, 0.05521),
    (0.04928, 0.06366),(0.06206, 0.07191),(0.07618, 0.07987),(0.09157, 0.08744),
    (0.10817, 0.09454),(0.12591, 0.10105),(0.14478, 0.10683),(0.16479, 0.11180),
    (0.18594, 0.11591),(0.20823, 0.11911),(0.23166, 0.12138),(0.25623, 0.12274),
    (0.28190, 0.12325),(0.30860, 0.12294),(0.33627, 0.12186),(0.36482, 0.12005),
    (0.39414, 0.11758),(0.42415, 0.11449),(0.45472, 0.11085),(0.48574, 0.10671),
    (0.51708, 0.10214),(0.54861, 0.09722),(0.58019, 0.09200),(0.61167, 0.08655),
    (0.64290, 0.08094),(0.67373, 0.07524),(0.70399, 0.06950),(0.73353, 0.06378),
    (0.76219, 0.05812),(0.78980, 0.05258),(0.81622, 0.04720),(0.84128, 0.04199),
    (0.86484, 0.03699),(0.88676, 0.03221),(0.90691, 0.02765),(0.92517, 0.02331),
    (0.94141, 0.01915),(0.95562, 0.01503),(0.96797, 0.01095),(0.97859, 0.00715),
    (0.98742, 0.00396),(0.99418, 0.00167),(0.99850, 0.00089),(1.00000, 0.00000),
]
_LOWER = [
    (0.00000, 0.00000),(0.00002,-0.00103),(0.00126,-0.00675),(0.00507,-0.01165),
    (0.01154,-0.01641),(0.02015,-0.02102),(0.03077,-0.02538),(0.04334,-0.02937),
    (0.05789,-0.03294),(0.07441,-0.03613),(0.09284,-0.03898),(0.11309,-0.04152),
    (0.13505,-0.04378),(0.15860,-0.04575),(0.18361,-0.04745),(0.20997,-0.04887),
    (0.23753,-0.05000),(0.26617,-0.05084),(0.29574,-0.05139),(0.32610,-0.05162),
    (0.35711,-0.05153),(0.38861,-0.05110),(0.42045,-0.05031),(0.45249,-0.04910),
    (0.48460,-0.04739),(0.51674,-0.04516),(0.54884,-0.04239),(0.58083,-0.03912),
    (0.61264,-0.03537),(0.64419,-0.03116),(0.67547,-0.02648),(0.70000,-0.02400),
    (0.72000,-0.02200),(0.74000,-0.02000),(0.76800,-0.01800),(0.79400,-0.01600),
    (0.82200,-0.01400),(0.85000,-0.01200),(0.88000,-0.01000),(0.89800,-0.00800),
    (0.92800,-0.00600),(0.95200,-0.00400),(0.98000,-0.00200),(1.00000, 0.00000),
]

upper = tuple((y * CHORD, z * CHORD) for y, z in _UPPER)
lower = tuple((y * CHORD, z * CHORD) for y, z in _LOWER)

def _interp_z(coords, y_tgt):
    for i in range(len(coords) - 1):
        y0, z0 = coords[i]; y1, z1 = coords[i+1]
        if y0 <= y_tgt <= y1:
            return z0 + (z1 - z0) * (y_tgt - y0) / (y1 - y0)
    return coords[-1][1]

ZS1T = _interp_z(upper, YSPAR1); ZS1B = _interp_z(lower, YSPAR1)
ZS2T = _interp_z(upper, YSPAR2); ZS2B = _interp_z(lower, YSPAR2)

# ── MODEL / MATERIAL / SECTIONS ───────────────────────────────────────────────
if MODEL_NAME in mdb.models: del mdb.models[MODEL_NAME]
mdb.Model(name=MODEL_NAME)
m = mdb.models[MODEL_NAME]

mat = m.Material(name='Aluminum')
mat.Elastic(table=((E, NU),))
mat.Density(table=((RHO,),))

for sec, t in (('SkinSec', SKIN_T), ('SparSec', SPAR_T), ('RibSec', RIB_T)):
    m.HomogeneousShellSection(name=sec, material='Aluminum',
        preIntegrate=OFF, thicknessType=UNIFORM, thickness=t,
        thicknessField='', nodalThicknessField='',
        idealization=NO_IDEALIZATION, integrationRule=SIMPSON, numIntPts=5)

# ── PART UTILITIES ────────────────────────────────────────────────────────────
def _assign(p, sec):
    p.SectionAssignment(region=regionToolset.Region(faces=p.faces[:]),
        sectionName=sec, offset=0., offsetType=MIDDLE_SURFACE,
        offsetField='', thicknessAssignment=FROM_SECTION)

def _mesh(p):
    p.seedPart(size=MESH_SIZE, deviationFactor=0.1, minSizeFactor=0.1)
    p.setElementType(regions=regionToolset.Region(faces=p.faces[:]),
        elemTypes=(mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD,
                                 hourglassControl=DEFAULT),
                   mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)))
    p.generateMesh()

def _sketch(p, plane, sz):
    dp = p.DatumPlaneByPrincipalPlane(principalPlane=plane, offset=0.)
    da = p.DatumAxisByPrincipalAxis(ZAXIS)
    tr = p.MakeSketchTransform(sketchPlane=p.datums[dp.id],
                                sketchUpEdge=p.datums[da.id],
                                sketchPlaneSide=SIDE1, origin=(0., 0., 0.))
    return m.ConstrainedSketch(name=p.name+'_sk', sheetSize=sz, transform=tr)

# ── PARTS ─────────────────────────────────────────────────────────────────────
def make_skin():
    # Closed profile: upper LE→TE, lower TE→LE — extruded along +X
    p  = m.Part(name='Skin', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    sk = _sketch(p, YZPLANE, CHORD * 3)
    sk.Spline(points=upper)
    sk.Spline(points=lower[::-1])
    p.BaseShellExtrude(sketch=sk, depth=SPAN)
    del m.sketches['Skin_sk']
    _assign(p, 'SkinSec')
    _mesh(p)
    return p

def make_spar(name, y_pos, zb, zt):
    # Rectangle in XZ plane (X=0→SPAN, Z=zb→zt), translated to Y=y_pos
    p  = m.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    sk = _sketch(p, XZPLANE, SPAN * 2)
    sk.rectangle(point1=(0., zb), point2=(SPAN, zt))
    p.BaseShell(sketch=sk)
    del m.sketches[name+'_sk']
    _assign(p, 'SparSec')
    _mesh(p)
    return p

def make_rib(name):
    # Inter-spar cross-section in YZ plane at X=0, translated in assembly
    zu_a = _interp_z(upper, YSPAR1); zl_a = _interp_z(lower, YSPAR1)
    zu_b = _interp_z(upper, YSPAR2); zl_b = _interp_z(lower, YSPAR2)
    u_tr = ([(YSPAR1, zu_a)] +
            [(y, z) for y, z in upper if YSPAR1 < y < YSPAR2] +
            [(YSPAR2, zu_b)])
    l_tr = ([(YSPAR1, zl_a)] +
            [(y, z) for y, z in lower if YSPAR1 < y < YSPAR2] +
            [(YSPAR2, zl_b)])
    p  = m.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    sk = _sketch(p, YZPLANE, CHORD * 2)
    sk.Spline(points=u_tr)
    sk.Spline(points=l_tr[::-1])
    sk.Line(point1=(YSPAR1, zu_a), point2=(YSPAR1, zl_a))
    sk.Line(point1=(YSPAR2, zu_b), point2=(YSPAR2, zl_b))
    p.BaseShell(sketch=sk)
    del m.sketches[name+'_sk']
    _assign(p, 'RibSec')
    _mesh(p)
    return p

skin_p  = make_skin()
fspar_p = make_spar('FrontSpar', YSPAR1, ZS1B, ZS1T)
rspar_p = make_spar('RearSpar',  YSPAR2, ZS2B, ZS2T)
rib_p   = make_rib('Rib')

# ── ASSEMBLY ──────────────────────────────────────────────────────────────────
asm = m.rootAssembly
asm.DatumCsysByDefault(CARTESIAN)

si  = asm.Instance(name='Skin-1',      part=skin_p,  dependent=ON)
fsi = asm.Instance(name='FrontSpar-1', part=fspar_p, dependent=ON)
rsi = asm.Instance(name='RearSpar-1',  part=rspar_p, dependent=ON)

# Spars created at Y=0, translate to chord station
asm.translate(instanceList=('FrontSpar-1',), vector=(SPAN, YSPAR1, 0.))
asm.translate(instanceList=('RearSpar-1',),  vector=(SPAN, YSPAR2, 0.))

rib_insts = []
for k, rx in enumerate(RIB_X):
    iname = 'Rib-%d' % (k + 1)
    ri = asm.Instance(name=iname, part=rib_p, dependent=ON)
    if rx > 0.:
        asm.translate(instanceList=(iname,), vector=(rx, 0., 0.))
    rib_insts.append((iname, rx, ri))

# ── STEP ──────────────────────────────────────────────────────────────────────
m.FrequencyStep(name='Frequency', previous='Initial',
                numEigen=N_EIGEN, eigensolver=LANCZOS,
                minEigen=None, maxEigen=None,
                blockSize=DEFAULT, maxBlocks=DEFAULT)

# ── TIE CONSTRAINTS ───────────────────────────────────────────────────────────
e = 0.05
zt_max = max(ZS1T, ZS2T); zb_min = min(ZS1B, ZS2B)
skin_faces = regionToolset.Region(faces=si.faces[:])

# Spar flange edges → skin faces
for spi, sy, zt, zb, tag in (
        (fsi, YSPAR1, ZS1T, ZS1B, 'FS'),
        (rsi, YSPAR2, ZS2T, ZS2B, 'RS')):
    m.Tie(name=tag+'_top', main=skin_faces,
          secondary=regionToolset.Region(edges=spi.edges.getByBoundingBox(
              xMin=-e, xMax=SPAN+e, yMin=sy-e, yMax=sy+e, zMin=zt-e, zMax=zt+e)),
          positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON)
    m.Tie(name=tag+'_bot', main=skin_faces,
          secondary=regionToolset.Region(edges=spi.edges.getByBoundingBox(
              xMin=-e, xMax=SPAN+e, yMin=sy-e, yMax=sy+e, zMin=zb-e, zMax=zb+e)),
          positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON)

# Rib boundary edges → skin faces + spar faces
for k, (iname, rx, ri) in enumerate(rib_insts):
    n = k + 1
    m.Tie(name='Rib%d_FS' % n,
          main=regionToolset.Region(faces=fsi.faces[:]),
          secondary=regionToolset.Region(edges=ri.edges.getByBoundingBox(
              xMin=rx-e, xMax=rx+e, yMin=YSPAR1-e, yMax=YSPAR1+e,
              zMin=ZS1B-0.1, zMax=ZS1T+0.1)),
          positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON)
    m.Tie(name='Rib%d_RS' % n,
          main=regionToolset.Region(faces=rsi.faces[:]),
          secondary=regionToolset.Region(edges=ri.edges.getByBoundingBox(
              xMin=rx-e, xMax=rx+e, yMin=YSPAR2-e, yMax=YSPAR2+e,
              zMin=ZS2B-0.1, zMax=ZS2T+0.1)),
          positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON)
    m.Tie(name='Rib%d_Skin' % n, main=skin_faces,
          secondary=regionToolset.Region(edges=ri.edges.getByBoundingBox(
              xMin=rx-e, xMax=rx+e, yMin=YSPAR1-0.1, yMax=YSPAR2+0.1,
              zMin=zb_min-0.1, zMax=zt_max+0.1)),
          positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON)

# ── BC — RP at root, kinematic coupling per part, encastre on RP ──────────────
root_rp = asm.ReferencePoint(point=(0., CHORD / 2., 0.))
rp_reg  = regionToolset.Region(
    referencePoints=(asm.referencePoints[root_rp.id],))

for inst, tag, y0, y1, z0, z1 in [
    (si,  'Skin', -0.1,        CHORD+0.1,   zb_min-0.1, zt_max+0.1),
    (fsi, 'FS',   YSPAR1-0.1,  YSPAR1+0.1,  ZS1B-0.1,   ZS1T+0.1),
    (rsi, 'RS',   YSPAR2-0.1,  YSPAR2+0.1,  ZS2B-0.1,   ZS2T+0.1),
]:
    edges = inst.edges.getByBoundingBox(
        xMin=-e, xMax=e, yMin=y0, yMax=y1, zMin=z0, zMax=z1)
    m.Coupling(name='Root_'+tag, controlPoint=rp_reg,
               surface=regionToolset.Region(edges=edges),
               influenceRadius=WHOLE_SURFACE, couplingType=KINEMATIC,
               localCsys=None, u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

m.EncastreBC(name='Cantilever', createStepName='Initial', region=rp_reg)
# ── OUTPUT + JOB ──────────────────────────────────────────────────────────────
m.fieldOutputRequests['F-Output-1'].setValues(variables=('U', 'UR',))

if JOB_NAME in mdb.jobs: del mdb.jobs[JOB_NAME]
mdb.Job(name=JOB_NAME, model=MODEL_NAME,
        description='Wing natural frequency',
        type=ANALYSIS, memory=90, memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
        echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
        scratch='', resultsFormat=ODB)

mdb.jobs[JOB_NAME].submit(consistencyChecking=OFF)
mdb.jobs[JOB_NAME].waitForCompletion()

print('Ready — %d ribs at X = %s' % (N_RIBS, RIB_X))
print('Submit: mdb.jobs["%s"].submit()' % JOB_NAME)