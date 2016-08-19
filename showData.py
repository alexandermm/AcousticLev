

import numpy as np
from mayavi import mlab



output = np.loadtxt("levData.txt")

N_levPoints = int(output[0, 0])
N_trans     = int(output[1, 0])

sf      = output[2, 0]
k       = output[3, 0]
P_0     = output[4, 0]
r       = output[5, 0]

K         = output[6, 0:2]
levPoints = output[7:(N_levPoints+7), :]

transPs = output[(N_levPoints+7):(N_levPoints+N_trans+7), :]
transNs = output[(N_levPoints+N_trans+7):(N_levPoints+2*N_trans+7), :]
phases  = output[(N_levPoints+N_trans*2+7):(N_levPoints+3*N_trans+7), 0]

print "\n"
print "N_levPoints = %d" % N_levPoints
print "N_trans     = %d" % N_trans

print "sf      = %f" % sf
print "k       = %f" % k
print "P_0     = %f" % P_0
print "r       = %f" % r
print "K:"
print K
print "levPoints:"
print levPoints
print "transPs shape:"
print transPs.shape
print "transNs shape:"
print transNs.shape
print "phases  shape:"
print phases.shape



#######################################################################################################################



#Transducer points with normal vectors
x = transPs[:,0]/sf
y = transPs[:,1]/sf 
z = transPs[:,2]/sf

u = transNs[:,0]/sf
v = transNs[:,1]/sf 
w = transNs[:,2]/sf  

vecs = mlab.pipeline.vector_scatter(x,y,z, u,v,w)
mlab.pipeline.vectors(vecs, scale_factor=20.0e-3)

transps = mlab.points3d(x,y,z, color=(0,0,1), scale_factor=4.0e-3)


#Phases
pu = u*phases
pv = v*phases
pw = w*phases
pha = mlab.pipeline.vector_scatter(x,y,z, pu,pv,pw)
mlab.pipeline.vectors(pha, scale_factor=20.0e-3)


#Levitation points
lx = levPoints[:,0]/sf
ly = levPoints[:,1]/sf 
lz = levPoints[:,2]/sf
levps = mlab.points3d(lx,ly,lz, color=(1,1,1), scale_factor=6.0e-3)


#Bounding box
gap    = 1.0e-2
height = 30.0e-2
p1 = np.array([np.min(x)-gap, np.max(y)+gap, -gap])
p2 = np.array([np.max(x)+gap, np.max(y)+gap, -gap])
p3 = np.array([np.max(x)+gap, np.min(y)-gap, -gap])
p4 = np.array([np.min(x)-gap, np.min(y)-gap, -gap])

p5 = np.array([p1[0],p1[1],height+gap])
p6 = np.array([p2[0],p2[1],height+gap])
p7 = np.array([p3[0],p3[1],height+gap])
p8 = np.array([p4[0],p4[1],height+gap])

bps = p1
bps = np.vstack([bps, p2])
bps = np.vstack([bps, p3])
bps = np.vstack([bps, p4])
bps = np.vstack([bps, p5])
bps = np.vstack([bps, p6])
bps = np.vstack([bps, p7])
bps = np.vstack([bps, p8])

x = bps[:,0]
y = bps[:,1] 
z = bps[:,2]
bounds = mlab.points3d(x,y,z, color=(1,1,1), scale_factor=5.0e-4)


#Show everything
mlab.axes()
mlab.outline()
mlab.show()

#######################################################################################################################




output = np.loadtxt("fieldData.txt")

numNodes = int(output.shape[0])/4

pts       = output[0:numNodes, :]
pressures = output[numNodes:numNodes*2,   0]
potential = output[numNodes*2:numNodes*3, 0]
forces    = output[numNodes*3:numNodes*4, :]

from tvtk.api import tvtk, write_data

x    = int(round((numNodes ** (1.0/3.0))+0.1))
dims = (x, x, x)
sg = tvtk.StructuredGrid(dimensions=dims, points=pts)
sg.point_data.scalars = (pressures.ravel())*sf
sg.point_data.scalars.name = 'Pressure magnitude'

sg.point_data.add_array((potential.ravel())/sf/sf)
sg.point_data.get_array(1).name = 'Gorkov potential'
sg.point_data.update()


#Get forces
sg.point_data.vectors = forces/sf
sg.point_data.vectors.name = 'Force field'


print "\nWriting .vtk file with pressure magnitude, Gorkov potential and force field\n"
write_data(sg, 'transducerField.vtk')




#######################################################################################################################




#Get path data
output = np.loadtxt("pathData.txt")

dt       = output[0,0]
pCenter  = output[1, 0:3]/sf
halfSide = output[2, 0:3]/sf
paths    = output[3:,  :]/sf

numParticles = paths.shape[0]/3
ncol         = paths.shape[1]


from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.gca(projection='3d')

plt.hold(True)

for n in range(numParticles):
	x = paths[3*n,:]
	y = paths[3*n+1,:]
	z = paths[3*n+2,:]
	ax.plot(x, y, z)


ax.set_xlabel('x-axis (m)')
ax.set_ylabel('y-axis (m)')
ax.set_zlabel('z-axis (m)')

plt.title('Particle paths')

axes = plt.gca()
axes.set_xlim([pCenter[0]-halfSide[0],pCenter[0]+halfSide[0]])
axes.set_ylim([pCenter[1]-halfSide[1],pCenter[1]+halfSide[1]])
axes.set_zlim([pCenter[2]-halfSide[2],pCenter[2]+halfSide[2]])


#Add lev points
ax.scatter(lx, ly, lz);

plt.show()


#######################################################################################################################


#Plot final points
output = np.loadtxt("positionData.txt")

diam_p       = output[0,0]
sf           = output[0,1]

numLevPoints   = int(output[0,2])
numFinalPoints = output.shape[0]-numLevPoints-3

pCenter  = output[1, :]/sf
halfSide = output[2, :]/sf

levPoints  = output[3:(numLevPoints+3), :]/sf
lastPoints = output[(numLevPoints+3):(output.shape[0]+1), :]/sf   


print "diam_p        = %f" % diam_p
print "sf            = %f" % sf
print "N_levPoints   = %d" % numLevPoints
print "N_finalPoints = %d" % numFinalPoints

print "pCenter:"
print pCenter
print "halfSide:"
print halfSide
print "levPoints:"
print levPoints
print "lastPoints:"
print lastPoints
print "\n"


fig = plt.figure()
ax = fig.gca(projection='3d')

plt.hold(True)

ax = fig.gca(projection='3d')


x = lastPoints[:, 0]
y = lastPoints[:, 1]
z = lastPoints[:, 2]
ax.scatter(x, y, z, c='red')


ax.set_xlabel('x-axis (m)')
ax.set_ylabel('y-axis (m)')
ax.set_zlabel('z-axis (m)')

plt.title('Final particle positions')

axes = plt.gca()
axes.set_xlim([pCenter[0]-halfSide[0],pCenter[0]+halfSide[0]])
axes.set_ylim([pCenter[1]-halfSide[1],pCenter[1]+halfSide[1]])
axes.set_zlim([pCenter[2]-halfSide[2],pCenter[2]+halfSide[2]])


#Add lev points
lx = levPoints[:, 0]
ly = levPoints[:, 1]
lz = levPoints[:, 2]
ax.scatter(lx, ly, lz, c='blue');

plt.show()



