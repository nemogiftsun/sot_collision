import dynamic_graph.sotcollision as sc
import numpy as np
from dynamic_graph.sot.core.matrix_util import matrixToTuple
ca = sc.SotCollision("sc")
ca.proximitySensorDistance.value = (0.07,)*8

p = np.zeros((7,8))
c = 0
for i in range(8):
    for j in range(7):
        c += 1
        p[j,i] = c
pose = matrixToTuple(p)


Mi = ((1, 0.0, 0, 0.1), 
      (0.0, 1.0, 0.0, 0.01),
      (0, 0.0, 0.3, 0.2),
      (0.0, 0.0, 0.0, 1.0))

Ji = ((1.0, 0.0, 0.0, 0.0, 0.089159, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
     (0.0, 1.0, 0.0, -0.089159, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
     (0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
     (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
     (0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
     (0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0))

ca.jointTransformation.value = Mi
ca.jointJacobian.value = Ji
ca.proximitySensorPose.value = pose
ca.setNumSkinSensors(8)
ca.collisionDistance.recompute(1)
ca.collisionJacobian.recompute(1)
      