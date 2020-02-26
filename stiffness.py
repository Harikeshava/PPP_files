import numpy as np 
import itertools as it
import math as mi 

def assemble_function(data,Assignment_matrixs,s,j):
    for name, age in data.items():    # for printing the value's key
        if age == s:
            i=name
            Assignment_matrixs[j][i]=1
            Assignment_matrixs[j+1][i+1]=1
    return(Assignment_matrixs)

def Material_routine(MU,lamda):
    #material routine will return the C_t matrixs
    C_elastic=np.array([[2*MU+lamda,lamda,0],[lamda,2*MU+lamda,0],[0,0,MU]])
    C_tangential=C_elastic
    return [C_tangential]

def Element_routine(Xe,Element_stiffness_matrixs,MU,lamda):
    #element routien will return the elements stiffness matrixs
    x_values=np.array([-0.57735,-0.57735,0.57735,-0.57735,0.57735,0.57735,-0.57735,0.57735])
     # C_elastic or C_tangential will be decided w.r.t the condition satisfied
    for i in range(0,8,2):
        x1=x_values[i]
        x2=x_values[i+1]
        derivative_x= np.array([[-(1-x2),(1-x2),(1+x2),-(1+x2)],[-(1-x1),-(1+x1),(1+x1),(1-x1)]]) * (1/4)
        Jacobin_matrixs = derivative_x @ Xe
        B_vector= np.linalg.inv(Jacobin_matrixs) @ derivative_x
        B=np.array([[B_vector[0][0],0,B_vector[0][1],0,B_vector[0][2],0,B_vector[0][3],0],[0,B_vector[1][0],0,B_vector[1][1],0,B_vector[1][2],0,B_vector[1][3]],[B_vector[1][0],B_vector[0][0],B_vector[1][1],B_vector[0][1],B_vector[1][2],B_vector[0][2],B_vector[1][3],B_vector[0][3]]])
        [C_tangential]=Material_routine(MU,lamda)
        Element_stiffness_matrixs= Element_stiffness_matrixs + (np.transpose(B) @ C_tangential @ B)* np.linalg.det(Jacobin_matrixs)
    return [Element_stiffness_matrixs]

#Elastic program for 2D Bilinear Element 
#Here, we are considering a point load as for building the basic strucutre for the required element
#internal parameters
yield_stress=0.0002
Youngs_modulus=0.2
Poissons_ratio=0.20
MU=(Youngs_modulus/(2*(1+Poissons_ratio)))
lamda=((Poissons_ratio*Youngs_modulus)/((1-2*Poissons_ratio)*(1+Poissons_ratio)))
L= eval(input('Enter the length of the plat\n'))
height_plate = eval(input('Enter the height of the plat\n'))
N=eval(input('Enter the number of elements in the x-direction\n'))
M=eval(input('Enter the number of elements in the y-direction\n'))
Le=L/N #element length
He=height_plate/M #element height
total_nodes= (N+1)*(M+1)
print(total_nodes)
Element_stiffness_matrixs=np.zeros((8,8))
Global_stiffness_matrixs=np.zeros((2*total_nodes,2*total_nodes))
print("Enter the forces value in newton for each node of interest\n")
print("Enter zero if no forces is applied on the node\n")
# The force part to be made into incremental wise after figuring the flow of the program
Xe=np.array([[0,0],[Le,0],[Le,He],[0,He]])
x_disp=np.array([[Le,0],[Le,0],[Le,0],[Le,0]])
y_disp=np.array([[0,He],[0,He],[0,He],[0,He]])
Node_Numbering= np.zeros(((M+1),(N+1)))
s=0
s1=0
s2=0
s3=0
s4=0
k=1
k1=0
data={}
for i in range(0,2*total_nodes,2):
    k1=k1+1
    data.update({i:k1})
for i in range(0,M+1):
    for j in range(0,N+1):
        Node_Numbering[i][j]=k 
        k=k+1
print(Node_Numbering)
d = np.sqrt((Le**2 + He**2))
Global_displacement=np.array([[0],[0],[Le*(mi.cos(30)-1)],[Le*mi.sin(30)],[d*mi.cos(30+45)-Le],[d*mi.sin(30+45)-He],[He*mi.cos(30+90)],[He*(mi.sin(30+90)-1)]])
for i in range(0,M):
    for j in range(0,N):
        Assignment_matrixs=np.zeros((8,2*total_nodes))
        s=i
        s1=Node_Numbering[s][j] #s1-first node of the local element
        s2=Node_Numbering[s][j+1] #s2-second node of the local element
        s3=Node_Numbering[s+1][j] #s3-third node of the local element
        s4=Node_Numbering[s+1][j+1] #s4-fourth node of the local element
        print('Element node:',s1,s2,s3,s4)
        Assignment_matrixs=assemble_function(data,Assignment_matrixs,s1,0)
        Assignment_matrixs=assemble_function(data,Assignment_matrixs,s2,2)
        Assignment_matrixs=assemble_function(data,Assignment_matrixs,s3,4)
        Assignment_matrixs=assemble_function(data,Assignment_matrixs,s4,6)
        if (j>0):
            Xe = Xe + x_disp
        [Element_stiffness_matrixs] = Element_routine(Xe,Element_stiffness_matrixs,MU,lamda)
        Global_stiffness_matrixs=Global_stiffness_matrixs + (np.transpose(Assignment_matrixs) @ Element_stiffness_matrixs @ Assignment_matrixs)
    Xe=Xe + y_disp       
print('The obtained stiffness matrixa is:\n')
print(Global_stiffness_matrixs)
F=Global_stiffness_matrixs @ Global_displacement
print(F)