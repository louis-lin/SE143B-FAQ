from abaqus import *
from abaqusConstants import *
from caeModules import *
import regionToolset
import mesh as meshMod

a, b, t   = 12.0, 5.0, 0.04
E, nu     = 10.0e6, 0.33
F_ref     = 1.0
seed_size = 0.20
tol       = 1e-5

Mdb()
mdl = mdb.models['Model-1']

sk = mdl.ConstrainedSketch(name='sk', sheetSize=20.0)
sk.rectangle(point1=(0.0, 0.0), point2=(a, b))
p = mdl.Part(name='Plate', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p.BaseShell(sketch=sk)

mdl.Material(name='Al').Elastic(table=((E, nu),))
mdl.HomogeneousShellSection(name='Sec', material='Al', thickness=t)
p.SectionAssignment(region=p.Set(name='All', faces=p.faces[:]), sectionName='Sec')
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
               numEigen=4, vectors=10, maxIterations=300)

# SSSS: u3=0 on all four edges; u1=0 at x=0 (no axial); u2=0 on RP (no rigid-body lateral)
mdl.DisplacementBC(name='BC_x0',   createStepName='Initial', region=Ex0, u1=SET, u3=SET)
mdl.DisplacementBC(name='BC_xa',   createStepName='Initial', region=Exa, u3=SET)
mdl.DisplacementBC(name='BC_y0',   createStepName='Initial', region=Ey0, u3=SET)
mdl.DisplacementBC(name='BC_yb',   createStepName='Initial', region=Eyb, u3=SET)
mdl.DisplacementBC(name='BC_RP',   createStepName='Initial', region=RP,  u2=SET)

mdl.ConcentratedForce(name='F', createStepName='Buckle', region=RP, cf1=-F_ref)

job = mdb.Job(name='PlateBuckle', model='Model-1', numCpus=2, numDomains=2)
job.writeInput()
job.submit()
job.waitForCompletion()
