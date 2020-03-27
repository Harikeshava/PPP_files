# Introduction:
The file contains information about the test cases used for the validation of the developed program.
The procedure to run the test case is mentioned in the user manual under the "Procedure to run the test cases" subsection.
# Test cases
## Rigid body Rotation:
Aim: The rigid rotation is a test done to verify that our system does not strain for a normal coordinate transformation.<br/>
Expected results:The system should be in equilibrium,$\sum_{i=1}^{N} F_{i}^{x}=0 and \sum_{i=1}^{N} F_{i}^{y}=0$<br/>
Obtained results: The obtained internal force matric is mentioned in figure-4.3.<br/>
Parameters used: The parameters used are mentioned below with the datas<br/>
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) --- The tension parameter is commented for this test case.<br/>
Global_displacement[2][0]= 4.909*10^-7 <br/>
Global_displacement[3][0]= 7.748*10^-8 <br/>
Global_displacement[4][0]= -4.909*10^-7 <br/>
Global_displacement[5][0]= -7.748*10**-8 <br/>
## Rigid body Rotation with tension:
Aim: The rigid rotation with tension, first the system is given a tension and then rotated. This again states that after rotation the system should not experience strain.<br/>
Expected results: the $F_{int}(U_{rig})=F_{int}(U_{rig+load\_applied}))$ <br/>
Obtained results: The obtained internal force matric is mentioned in figure- 4.4 and 4.5. <br/>
Parameters used: <br/>
tension_disp=np.array([[0,0],[4.909*10^-9,7.751*10^-10],[0,0],[4.909*10^-9,-7.751*10^-10]])
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) + tension_disp
Global_displacement[2][0]= 4.909*10^-7 
Global_displacement[3][0]= 7.748*10^-8
Global_displacement[4][0]= -4.909*10^-7
Global_displacement[5][0]= -7.748*10^-8
## Rigid body Rotation with compression:
Aim: The rigid rotation with compression, first the system is given a tension and then rotated. This again states that after rotation the system should not experience strain. <br/>
Expected results: the $F_{int}(U_{rig})=F_{int}(U_{rig+load\_applied}))$ <br/>
Obtained results: The obtained internal force is mentioned in figure- 4.6 and 4.7. <br/>
Parameters used: <br/>
tension_disp=np.array([[0,0],[-4.909*10^-9,-7.751*10^-10],[0,0],[-4.909*10^-9,7.751*10^-10]]) <br/>
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) + tension_disp <br/>
Global_displacement[2][0]= 4.909*10^-7 <br/>
Global_displacement[3][0]= 7.748*10^-8 <br/>
Global_displacement[4][0]= -4.909*10^-7 <br/>
Global_displacement[5][0]= -7.748*10^-8 <br/>
## Rigid body Rotation with shear:
Aim: The rigid rotation with shear, first the system is given a tension and then rotated. This again states that after rotation the system should not experience strain. <br/>
Expected results: the $F_{int}(U_{rig})=F_{int}(U_{rig+load\_applied}))$ <br/>
Obtained results: The obtained internal force is mentioned in figure- 4.8 and 4.9. <br/>
Parameters used: <br/>
tension_disp=np.array([[0,0],[5.55*10^-9,1.88*10^-10],[0,0],[-6.33*10^-9,1.97*10^-10]]) <br/>
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) + tension_disp <br/>
Global_displacement[2][0]= 4.909*10^-7  <br/>
Global_displacement[3][0]= 7.748*10^-8 <br/>
Global_displacement[4][0]= -4.909*10^-7 <br/>
Global_displacement[5][0]= -7.748*10^-8 <br/>
