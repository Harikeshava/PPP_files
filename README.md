# Introduction
The file contains information about the test cases used for the validation of the developed program.
The procedure to run the test case is mentioned in the user manual under the "Procedure to run the test cases" subsection.
# Macro Test cases
## Rigid body Rotation:
Aim: The rigid rotation is a test done to verify that our system does not strain for a normal coordinate transformation.<br/>
Expected results:The system should be in equilibrium,that is the forces in x and y direction should be equal to zero as mentioned in the report.<br/>
Obtained results: The obtained internal force matric is mentioned in figure-4.3.<br/>
Parameters used: The parameters used are mentioned below with the datas. The parameters Le and He are assigned in the code using the input parameters.<br/>
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) --- The tension parameter is commented for this test case.<br/>
Global_displacement[2][0]= 4.909*10^-7 <br/>
Global_displacement[3][0]= 7.748*10^-8 <br/>
Global_displacement[4][0]= -4.909*10^-7 <br/>
Global_displacement[5][0]= -7.748*10^-8 <br/>
## Rigid body Rotation with tension:
Aim: The rigid rotation with tension, first the system is given a tension and then rotated. This again states that after rotation the system should not experience strain.<br/>
Expected results: the F<sub>int</sub>(U<sub>rig</sub>)=F<sub>int</sub>(U<sub>rig+load\_applied</sub>))  <br/>
Obtained results: The obtained internal force matric is mentioned in figure- 4.4 and 4.5. <br/>
Parameters used: <br/>
tension_disp=np.array([[4.909*10^-9,0],[4.909*10^-9,0],[4.909*10^-9,0],[4.909*10^-9,0]]) <br/>
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) + tension_disp <br/>
Global_displacement[2][0]= 4.909×10^-7 <br/>
Global_displacement[3][0]= 7.748×10^-8 <br/>
Global_displacement[4][0]= -4.909×10^-7 <br/>
Global_displacement[5][0]= -7.748×10^-8 <br/>
## Rigid body Rotation with compression:
Aim: The rigid rotation with compression, first the system is given a tension and then rotated. This again states that after rotation the system should not experience strain. <br/>
Expected results: the F<sub>int</sub>(U<sub>rig</sub>)=F<sub>int</sub>(U<sub>rig+load\_applied</sub>))  <br/>
Obtained results: The obtained internal force is mentioned in figure- 4.6 and 4.7. <br/>
Parameters used: <br/>
tension_disp=np.array([[-4.909*10^-9,0],[-4.909*10^-9,0],[-4.909*10^-9,0],[-4.909*10^-9,0]]) <br/>
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) + tension_disp <br/>
Global_displacement[2][0]= 4.909*10^-7 <br/>
Global_displacement[3][0]= 7.748*10^-8 <br/>
Global_displacement[4][0]= -4.909*10^-7 <br/>
Global_displacement[5][0]= -7.748*10^-8 <br/>
## Rigid body Rotation with shear:
Aim: The rigid rotation with shear, first the system is given a tension and then rotated. This again states that after rotation the system should not experience strain. <br/>
Expected results: the F<sub>int</sub>(U<sub>rig</sub>)=F<sub>int</sub>(U<sub>rig+load\_applied</sub>))   <br/>
Obtained results: The obtained internal force is mentioned in figure- 4.8 and 4.9. <br/>
Parameters used: <br/>
tension_disp= np.array([[4.909*10**-7,0],[-3.27266*10**-8,0],[3.27266*10**-8,0],[-4.909*10**-7,0]]) <br/>
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) + tension_disp <br/>
Global_displacement[2][0]= 4.909×10^-7  <br/>
Global_displacement[3][0]= 7.748×10^-8 <br/>
Global_displacement[4][0]= -4.909×10^-7 <br/>
Global_displacement[5][0]= -7.748×10^-8 <br/>
#Micro Test case: 
The macro strain parameters for all the cases is as displayed below.<br/>
Strain_macro=np.array([[5.16754*10**-7],[-1.0335*10**-7],[-4.268*10**-7]])<br/>