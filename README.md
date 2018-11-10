# BEM_flow_simulaiton
Comutational fluid mechanics on MATLAB :

The goal is to simulate any kind of flow near or around any kind of obstacle using the Boundary Elements Methods

1. Potential flow around cylindric obstacle &amp; near uniform war using BEM 

The hypotheses of non viscous potential flow are applyed to the navier stokes equations (Euler equations), using the boundary elements methods, this code simulates the flow around a dimensionally stable cylindric obstacle that is near a uniform vertical wall. The system is free from any boundaries on the other directions. This code was the subject of a master 1 internship in fluid mechanics, and as a first brick of a big project, which is to simulate viscous flows around complexe non uniform obstacles, it is free to be continued.

The class "pot_flow_class" simulates this potential flow around the cylindric 2D obstacle (dimensionless radius r=1) near the unifrom vertical wall (distance between the two obstacle is H). The flow velocity is U=1 (dimensionless). The second layer of the boundary elements methods is used. A code of the exact solution of such flow allows to verify the no vertical wall simulation. 
