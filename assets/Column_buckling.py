# -*- coding: utf-8 -*-
"""
Linear buckling: 5x0.04 in Al strip, pin-pin, L=12 in
Hand: Pcr = pi^2*E*I/L^2 ~ 18.3 lb  (I = 5*0.04^3/12 = 2.67e-5 in^4)

FIX: Reference load must live in the BuckleStep itself, not a prior StaticStep.
     StaticStep preloads are only needed for physical preloads (gravity, etc.)
     that alter the base stress state before buckling is assessed.

Run: abaqus cae noGUI=strip_buckling.py
"""
from abaqus import *
from abaqusConstants import *
import regionToolset, mesh, math

NAME        = 'StripBuckling'
E, NU       = 10.0e6, 0.33
L, b_w, t   = 12.0, 5.0, 0.04
N_EIGEN     = 3

mdb.Model(name=NAME)
m = mdb.models[NAME]

if 'Model-1' in mdb.models:
    del mdb.models['Model-1']

mat = m.Material(name='Al')
mat.Elastic(table=((E, NU),))

# a=t along n1=Y, b=b_w along n2=Z -> I_11 = b_w*t^3/12 (weak axis, min Pcr)
m.RectangularProfile(name='Prof', a=t, b=b_w)
m.BeamSection(name='Sec', profile='Prof', material='Al',
              poissonRatio=0.0, integration=DURING_ANALYSIS, beamShape=CONSTANT)

p = m.Part(name='Strip', dimensionality=THREE_D, type=DEFORMABLE_BODY)
p.WirePolyLine(points=((0.,0.,0.),(L,0.,0.)), mergeType=IMPRINT, meshable=ON)
eReg = regionToolset.Region(edges=p.edges[:])
p.SectionAssignment(region=eReg, sectionName='Sec')
p.assignBeamSectionOrientation(region=eReg, method=N1_COSINES, n1=(0.,1.,0.))

p.seedPart(size=0.25)
p.setElementType(regions=eReg,
                 elemTypes=(mesh.ElemType(elemCode=B31, elemLibrary=STANDARD),))
p.generateMesh()

a = m.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
inst = a.Instance(name='Strip-1', part=p, dependent=ON)

# Buckle step directly after Initial — no StaticStep needed
m.BuckleStep(name='Buckle', previous='Initial',
             numEigen=N_EIGEN, vectors=2*N_EIGEN+1, maxIterations=30)

r0 = a.Set(name='Left',  vertices=inst.vertices.findAt(((0.,0.,0.),)))
rL = a.Set(name='Right', vertices=inst.vertices.findAt(((L, 0.,0.),)))

# Left pin:  U1=U2=U3=0, rotations free
m.DisplacementBC(name='BC_L', createStepName='Initial', region=r0,
                 u1=0.0, u2=0.0, u3=0.0)
# Right pin: U2=U3=0, U1 free (load end)
m.DisplacementBC(name='BC_R', createStepName='Initial', region=rL,
                 u2=0.0, u3=0.0)

# 1 lb reference load IN BuckleStep -> eigenvalue = Pcr in lb
m.ConcentratedForce(name='Pref', createStepName='Buckle', region=rL, cf1=-1.0)

mdb.Job(name=NAME, model=NAME, numCpus=1, memory=90,
        memoryUnits=PERCENTAGE, resultsFormat=ODB)
mdb.jobs[NAME].submit(consistencyChecking=OFF)
mdb.jobs[NAME].waitForCompletion()

# Open ODB and render beam cross-sections
odb = session.openOdb(name=NAME + '.odb')
vp  = session.viewports[session.currentViewportName]
vp.setValues(displayedObject=odb)
vp.odbDisplay.display.setValues(plotState=(DEFORMED,))
vp.odbDisplay.basicOptions.setValues(renderBeamProfiles=ON, beamScaleFactor=1)
vp.odbDisplay.setPrimaryVariable(variableLabel='U', outputPosition=NODAL,
                                  refinement=(INVARIANT, 'Magnitude'))
