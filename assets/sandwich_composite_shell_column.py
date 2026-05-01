from abaqus import *
from abaqusConstants import *
from section import SectionLayer
import regionToolset, mesh
import math

E_AL, NU_AL = 10.0e6, 0.33
E_FM, NU_FM = 10.0E6, 0.3    # G = E/2(1+nu) ~ 1154 psi  ← shear-capable
b_w, L, t_s = 5.0, 12.0, 0.04
MM      = 1.0 / 25.4
EPS     = 1e-4
N_EIGEN = 5

mdb.Model(name='_')
if 'Model-1' in mdb.models:
    del mdb.models['Model-1']

for i, t_c_mm in enumerate((2,3,4)):
    t_c  = t_c_mm * MM
    NAME = 'Sandwich_%dmm_shell' % t_c_mm

    mdb.Model(name=NAME)
    if i == 0 and '_' in mdb.models:
        del mdb.models['_']
    m = mdb.models[NAME]

    m.Material(name='Al').Elastic(table=((E_AL, NU_AL),))
    m.Material(name='Foam').Elastic(table=((E_FM, NU_FM),))

    m.CompositeShellSection(
        name='Sec', preIntegrate=OFF, symmetric=False,
        idealization=NO_IDEALIZATION, thicknessType=UNIFORM,
        poissonDefinition=DEFAULT, temperature=GRADIENT, useDensity=OFF,
        layup=(
            SectionLayer(thickness=t_s/2, material='Al',   orientAngle=0., numIntPts=3),
            SectionLayer(thickness=t_c, material='Foam', orientAngle=0., numIntPts=1),
            SectionLayer(thickness=t_s/2, material='Al',   orientAngle=0., numIntPts=3),
        )
    )

    sk = m.ConstrainedSketch(name='sk', sheetSize=20.)
    sk.rectangle(point1=(0., 0.), point2=(b_w, L))
    p = m.Part(name='Panel', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseShell(sketch=sk)

    fReg = regionToolset.Region(faces=p.faces[:])
    p.SectionAssignment(region=fReg, sectionName='Sec', offset=0.,
                        offsetType=MIDDLE_SURFACE, offsetField='',
                        thicknessAssignment=FROM_SECTION)
    p.seedPart(size=0.25)
    p.setElementType(regions=fReg,
                     elemTypes=(mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD),))
    p.generateMesh()

    asm  = m.rootAssembly
    asm.DatumCsysByDefault(CARTESIAN)
    inst = asm.Instance(name='Panel-1', part=p, dependent=ON)

    # vectors must exceed numEigen — Lanczos basis must span the eigenspace
    m.BuckleStep(name='Buckle', previous='Initial',
                 numEigen=N_EIGEN, vectors=2*N_EIGEN+1, maxIterations=10000)

    e0 = inst.edges.getByBoundingBox(yMin=-EPS, yMax=EPS)
    eL = inst.edges.getByBoundingBox(yMin=L-EPS, yMax=L+EPS)

    rp0_obj = asm.ReferencePoint(point=(b_w/2., 0., 0.))
    rpL_obj = asm.ReferencePoint(point=(b_w/2., L,  0.))
    rp0 = asm.referencePoints[rp0_obj.id]
    rpL = asm.referencePoints[rpL_obj.id]

    set_rp0 = asm.Set(name='RP_Bottom', referencePoints=(rp0,))
    set_rpL = asm.Set(name='RP_Top',    referencePoints=(rpL,))
    surf_e0 = asm.Surface(name='Surf_E0', side1Edges=e0)
    surf_eL = asm.Surface(name='Surf_EL', side1Edges=eL)

    # DISTRIBUTING on both ends — weighted kinematic, avoids rigid endplate
    # artifact that KINEMATIC causes (stress concentration + spurious local modes)
    for coup_name, ctrl, surf in (
        ('Coup_0', set_rp0, surf_e0),
        ('Coup_L', set_rpL, surf_eL),
    ):
        m.Coupling(name=coup_name, controlPoint=ctrl, surface=surf,
                   couplingType=DISTRIBUTING, influenceRadius=WHOLE_SURFACE,
                   u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)

    # Bottom RP: pin (u1=u2=u3=0, ur2=ur3=0, ur1 free → column rotation free)
    m.DisplacementBC(name='BC_Bottom', createStepName='Initial', region=set_rp0,
                     u1=0., u2=0., u3=0., ur2=0., ur3=0.)

    # Top RP: pin with axial freedom (u2 free to compress)
    m.DisplacementBC(name='BC_Top', createStepName='Initial', region=set_rpL,
                     u1=0., u3=0., ur2=0., ur3=0.)

    # Total reference load = 1 lb (eigenvalue = Pcr in lb)
    m.ConcentratedForce(name='Load_P', createStepName='Buckle',
                        region=set_rpL, cf2=-1.0)

    mdb.Job(name=NAME, model=NAME, numCpus=1, memory=90,
            memoryUnits=PERCENTAGE, resultsFormat=ODB)
    mdb.jobs[NAME].submit(consistencyChecking=OFF)
    mdb.jobs[NAME].waitForCompletion()

    # Print theoretical I and Pcr for each core (add inside the loop, after t_c is defined)

    t_skin = t_s / 2.0
    d      = (t_c + t_skin) / 2.0
    I11    = 2 * b_w * (t_skin**3 / 12.0 + t_skin * d**2)
    Pcr_hand = math.pi**2 * E_AL * I11 / L**2
    print('Core %dmm | I11=%.4e in^4 | Hand Pcr=%.1f lb' % (t_c_mm, I11, Pcr_hand))