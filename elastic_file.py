import numpy as np 
import itertools as it
def u_arrangement(data,Global_displacement,s,u_element,t):
    for name, age in data.items():    # for printing the value's key
        if age == s:
            i=name
            u_element[t][0]=Global_displacement[i][0]
            u_element[t+1][0]=Global_displacement[i+1][0]
    return u_element

def assemble_function(data,Assignment_matrixs,s,j):
    for name, age in data.items():    # for printing the value's key
        if age == s:
            i=name
            Assignment_matrixs[j][i]=1
            Assignment_matrixs[j+1][i+1]=1
    return(Assignment_matrixs)

def Material_routine(MU,lamda,u_element,B):
    #material routine will return the C_t matrixs
    C_elastic=np.array([[2*MU+lamda,lamda,0],[lamda,2*MU+lamda,0],[0,0,MU]])
    C_tangential=C_elastic 
    Strain_element=np.zeros((3,1))
    Strain_element= B @ u_element 
    print("Strain_element:\n",Strain_element)
    stress_element= C_tangential @ Strain_element
    return [C_tangential,stress_element]

def Element_routine(Xe,Element_stiffness_matrixs,MU,lamda,u_element,F_int_Element,Le,thickness_plate):
    #element routien will return the elements stiffness matrixs
    x_values=np.array([-0.57735,-0.57735,0.57735,-0.57735,-0.57735,0.57735,0.57735,0.57735])
    Element_stiffness_matrixs=np.zeros((8,8))
    F_int_Element=np.zeros((8,1))
    for i in range(0,8,2):
        x1=x_values[i]
        x2=x_values[i+1]
        derivative_x= np.array([[-(1-x2),(1-x2),-(1+x2),(1+x2)],[-(1-x1),-(1+x1),(1-x1),(1+x1)]]) * (1/4)
        Jacobin_matrixs = derivative_x @ Xe
        B_vector= np.linalg.inv(Jacobin_matrixs) @ derivative_x
        B=np.array([[B_vector[0][0],0,B_vector[0][1],0,B_vector[0][2],0,B_vector[0][3],0],[0,B_vector[1][0],0,B_vector[1][1],0,B_vector[1][2],0,B_vector[1][3]],[B_vector[1][0],B_vector[0][0],B_vector[1][1],B_vector[0][1],B_vector[1][2],B_vector[0][2],B_vector[1][3],B_vector[0][3]]])
        [C_tangential,stress_element]=Material_routine(MU,lamda,u_element,B) 
        Element_stiffness_matrixs= Element_stiffness_matrixs + (np.transpose(B) @ C_tangential @ B)* np.linalg.det(Jacobin_matrixs)*thickness_plate
        F_int_Element= F_int_Element + (np.transpose(B) @ stress_element ) * np.linalg.det(Jacobin_matrixs)*thickness_plate
    F_int_Element_1= Element_stiffness_matrixs @ u_element
    #if(np.linalg.norm(F_int_Element)==np.linalg.norm(F_int_Element_1)):
        #print("yes")
    return [Element_stiffness_matrixs,F_int_Element]

#Elastic program for 2D Bilinear Element 
#Here, we are considering a point load as for building the basic strucutre for the required element
#internal parameters
yield_stress=70*10**-6
Youngs_modulus=210E9 #N/meters
Poissons_ratio=0.30
MU=(Youngs_modulus/(2*(1+Poissons_ratio)))
lamda=((Poissons_ratio*Youngs_modulus)/((1-2*Poissons_ratio)*(1+Poissons_ratio)))
print(MU)
print(lamda)
L= eval(input('Enter the length of the plat in meters\n'))
height_plate = eval(input('Enter the height of the plat in meters\n'))
thickness_plate = eval(input("Enter the thickness of the plate meters\n"))
N=eval(input('Enter the number of elements in the x-direction\n'))
M=eval(input('Enter the number of elements in the y-direction\n'))
Le=L/N #element length
He=height_plate/M #element height
total_nodes= (N+1)*(M+1)
print(total_nodes)
delta_u=np.zeros((2*total_nodes,1))
Element_stiffness_matrixs=np.zeros((8,8))
Global_F_external=np.zeros((2*total_nodes,1))
Global_plastic_strain=np.zeros((3,M*N))
Global_displacement=np.zeros((2*total_nodes,1))
u_element=np.zeros((8,1))
F_int_Element=np.zeros((8,1))
#print("Enter the forces value in newton for each node of interest\n")
#print("Enter zero if no forces is applied on the node\n")
# The force part to be made into incremental wise after figuring the flow of the program
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]])
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
R_delta_u=np.ones((2*total_nodes,1))
for i in range(0,2*total_nodes,2):
    k1=k1+1
    data.update({i:k1})
for i in range(0,M+1):
    for j in range(0,N+1):
        Node_Numbering[i][j]=k 
        k=k+1
print(Node_Numbering)
count=0
while(np.linalg.norm(R_delta_u,np.inf) > (0.005*np.linalg.norm(Global_displacement,np.inf))):
    count=count+1
    Global_stiffness_matrixs=np.zeros((2*total_nodes,2*total_nodes))
    Global_F_internal=np.zeros((2*total_nodes,1))
    print(np.linalg.norm(R_delta_u,np.inf))
    print(np.linalg.norm(Global_displacement,np.inf))
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
            u_element=u_arrangement(data,Global_displacement,s1,u_element,0)
            u_element=u_arrangement(data,Global_displacement,s2,u_element,2)
            u_element=u_arrangement(data,Global_displacement,s3,u_element,4)
            u_element=u_arrangement(data,Global_displacement,s4,u_element,6)
            if (j>0):
                Xe = Xe + x_disp
            [Element_stiffness_matrixs,F_int_Element] = Element_routine(Xe,Element_stiffness_matrixs,MU,lamda,u_element,F_int_Element,Le,thickness_plate)
            Global_stiffness_matrixs=Global_stiffness_matrixs + (np.transpose(Assignment_matrixs) @ Element_stiffness_matrixs @ Assignment_matrixs)
            F=(np.transpose(Assignment_matrixs) @ F_int_Element)
            Global_F_internal=Global_F_internal + F
        if((M*N)>1):
            Xe=Xe + y_disp       
    print('The obtained stiffness matrixa is:\n')
    print(Global_stiffness_matrixs)
    print("Internal_force:\n")
    print(Global_F_internal)
    Global_F_external[2][0]=1000
    Global_F_external[6][0]=1000
    G=Global_F_internal-Global_F_external
    print("G:",G)
    Reduced_Global_stiffness_matrix=Global_stiffness_matrixs
    Reduced_displacement=Global_displacement
    Reduced_G=G
    A=[]
    #Reduction of matirxs sizes
    for j in range(0,M+1,1):    
        for name, age in data.items():    # for printing the value's key
            if age == Node_Numbering[j][0]:
                i=name
                A.append(i)
                A.append(i+1)
    A=np.asarray(A)
    print(A)
    if(count==1):
        R_delta_u=np.delete(delta_u,A,axis=0)
    Reduced_Global_stiffness_matrix= np.delete(Reduced_Global_stiffness_matrix,A,axis=0)
    Reduced_Global_stiffness_matrix= np.delete(Reduced_Global_stiffness_matrix,A,axis=1)
    Reduced_G=np.delete(Reduced_G,A,axis=0)
    Reduced_displacement=np.delete(Reduced_displacement,A,axis=0)
    #Reduced_F_ext=np.delete(Global_F_external,A,axis=0)
    #calculation part
    print("reduced k matrixs",Reduced_Global_stiffness_matrix)
    K_inv=np.linalg.inv(Reduced_Global_stiffness_matrix)
    #R_delta_u= K_inv @ Reduced_F_ext
    #print(R_delta_u)
    print("Reduced G",Reduced_G)
    R_delta_u=(np.linalg.inv(Reduced_Global_stiffness_matrix)) @ Reduced_G
    print("r_DELTA_U",R_delta_u)
    Reduced_displacement=Reduced_displacement - R_delta_u
    for i in range(0,len(A),1):
        Reduced_displacement=np.insert(Reduced_displacement,A[i],Global_displacement[A[i]])
        #R_delta_u=np.insert(R_delta_u,A[i],delta_u[i])
    Global_displacement=(Reduced_displacement.reshape(2*total_nodes,1))
print("displacement:\n",Global_displacement)
print("Iteration number:",count)