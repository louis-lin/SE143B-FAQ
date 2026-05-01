from abaqus import *
from abaqusConstants import *
from caeModules import *
from section import SectionLayer
import regionToolset
import mesh as meshMod

a, b, t_s = 12.0, 5.0, 0.04
t_c_mm    = 3
t_c       = t_c_mm / 25.4
E_AL, NU_AL = 10.0e6, 0.33
E_FM, NU_FM = 3000.0, 0.3
F_ref     = 1.0
seed_size = 0.20
tol       = 1e-5
NAME      = 'SandwichPlate_3mm'

Mdb()
mdb.Model(name=NAME)
if 'Model-1' in mdb.models:
    del mdb.models['Model-1']
mdl = mdb.models[NAME]

sk = mdl.ConstrainedSketch(name='sk', sheetSize=20.0)
sk.rectangle(point1=(0.0, 0.0), point2=(a, b))
p = mdl.Part(name='Plate', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p.BaseShell(sketch=sk)

mdl.Material(name='Al').Elastic(table=((E_AL, NU_AL),))
mdl.Material(name='Foam').Elastic(table=((E_FM, NU_FM),))
mdl.CompositeShellSection(
    name='Sec', preIntegrate=OFF, symmetric=False,
    idealization=NO_IDEALIZATION, thicknessType=UNIFORM,
    poissonDefinition=DEFAULT, temperature=GRADIENT, useDensity=OFF,
    layup=(
        SectionLayer(thickness=t_s/2, material='Al',   orientAngle=0., numIntPts=3),
        SectionLayer(thickness=t_c,   material='Foam', orientAngle=0., numIntPts=1),
        SectionLayer(thickness=t_s/2, material='Al',   orientAngle=0., numIntPts=3),
    )
)

fReg = regionToolset.Region(faces=p.faces[:])
p.SectionAssignment(region=fReg, sectionName='Sec', offset=0.,
                    offsetType=MIDDLE_SURFACE, offsetField='',
                    thicknessAssignment=FROM_SECTION)
p.setElementType(regions=(p.faces,),
    elemTypes=(meshMod.ElemType(elemCode=S4R, elemLibrary=STANDARD),))
p.seedPart(size=seed_size)
p.generateMesh()

asm  = mdl.rootAssembly
inst = asm.Instance(name='I', part=p, dependent=ON)

def edge_set(name, xmin, ymin, xmax, ymax):
    e = inst.edges.getByBoundingBox(xmin-tol, ymin-tol, -tol,
                                    xmax+tol, ymax+tol,  tol)
    return asm.Set(name=name, edges=e)

Ex0 = edge_set('Ex0', 0., 0., 0., b)
Exa = edge_set('Exa', a,  0., a,  b)
Ey0 = edge_set('Ey0', 0., 0., a,  0.)
Eyb = edge_set('Eyb', 0., b,  a,  b)
asm.Surface(name='Sxa', side1Edges=Exa.edges)

rp_id = asm.ReferencePoint(point=(a, b/2.0, 0.0)).id
RP    = asm.Set(name='RP', referencePoints=(asm.referencePoints[rp_id],))

mdl.Coupling(name='cpl', controlPoint=RP, surface=asm.surfaces['Sxa'],
             influenceRadius=WHOLE_SURFACE, couplingType=DISTRIBUTING,
             weightingMethod=UNIFORM,
             u1=ON, u2=ON, u3=ON, ur1=OFF, ur2=OFF, ur3=OFF)

mdl.BuckleStep(name='Buckle', previous='Initial',
               numEigen=4, vectors=10, maxIterations=1000)

mdl.DisplacementBC(name='BC_x0', createStepName='Initial', region=Ex0, u1=SET, u3=SET)
mdl.DisplacementBC(name='BC_xa', createStepName='Initial', region=Exa, u3=SET)
mdl.DisplacementBC(name='BC_y0', createStepName='Initial', region=Ey0, u3=SET)
mdl.DisplacementBC(name='BC_yb', createStepName='Initial', region=Eyb, u3=SET)
mdl.DisplacementBC(name='BC_RP', createStepName='Initial', region=RP,  u2=SET)

mdl.ConcentratedForce(name='F', createStepName='Buckle', region=RP, cf1=-F_ref)

job = mdb.Job(name=NAME, model=NAME, numCpus=2, numDomains=2)
job.writeInput()
job.submit()
job.waitForCompletion()
