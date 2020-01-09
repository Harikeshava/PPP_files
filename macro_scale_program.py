import numpy as np 
import itertools as it
def u_arrangement(data,Global_displacement,s,u_element,t): # The function for arranging  the displacement element accessed
    for name, age in data.items():    # for printing the value's key
        if age == s:
            i=name
            u_element[t][0]=Global_displacement[i][0]
            u_element[t+1][0]=Global_displacement[i+1][0]
    return u_element

def assemble_function(data,Assignment_matrixs,s,j): # The function to assemble the assignment matrixs
    for name, age in data.items():    # for printing the value's key
        if age == s:
            i=name
            Assignment_matrixs[j][i]=1
            Assignment_matrixs[j+1][i+1]=1
    return(Assignment_matrixs)

def Material_routine(MU,lamda,u_element,B,h,B_ps,C_al,P,j,yield_stress): # Material routine 
    #material routine will return the C_t matrixs,stress at gauss points and the internal state variables
    Norm_vector=np.zeros((6,1))
    P_sym=np.array([[(2/3),(-1/3),(-1/3),0,0,0],[(-1/3),(2/3),(-1/3),0,0,0],[(-1/3),(-1/3),(2/3),0,0,0],[0,0,0,(1/2),0,0],[0,0,0,0,(1/2),0],[0,0,0,0,0,(1/2)]])
    Strain_element=np.zeros((3,1))
    Strain_element= B @ u_element
    for i in range(2,5,1):
        Strain_element= np.insert(Strain_element,i,0,axis=0)
    #print("Strain_element:\n",Strain_element)
    trace_strain= Strain_element[0][0]+Strain_element[1][0]+Strain_element[2][0]
    Identity_matrixs=np.array([[1],[1],[1],[0],[0],[0]])
    deviatoric_strain= Strain_element - (trace_strain*Identity_matrixs)/3
    for i in range(3,6,1):
        deviatoric_strain[i][0]=(deviatoric_strain[i][0])/2
    plastic_strain=(B_ps[j].reshape(1,6))
    #print(plastic_strain)
    alpha=C_al[j][0]
    deviatoric_stress= 2*MU*(deviatoric_strain-np.transpose(plastic_strain))
    deviatoric_stress_trial=deviatoric_stress
    s=0
    s1=1
    for i in range(0,6,1):
        if(i>2):
            s1=2
        s=s+(deviatoric_stress_trial[i][0]*deviatoric_stress_trial[i][0])*s1
    mag_deviatoric=(s)**(1/2)
    if(mag_deviatoric > 0):
        Norm_vector=deviatoric_stress/mag_deviatoric
    drag_stress=h*alpha
    plastic_fn=mag_deviatoric-(2/3)**(1/2)*(yield_stress + drag_stress)
    if(plastic_fn<10**-6):
        print("Elastic calculation")
        deviatoric_stress=deviatoric_stress_trial
        deviatoric_tangent_stiffness= 2*MU*P_sym
        plastic_strain=np.transpose(plastic_strain)
        alpha=np.transpose(alpha)
    else:
        print("Plastic calculation")
        plastic_corrector= (plastic_fn)/(2*MU + (2/3)*h)
        deviatoric_stress = deviatoric_stress_trial - 2*MU*Norm_vector*plastic_corrector
        plastic_strain= np.transpose(plastic_strain) + plastic_corrector*Norm_vector
        alpha = np.transpose(alpha) + (2/3)**(1/2)*(plastic_corrector)
        S2 = (plastic_fn/mag_deviatoric)
        S3 = (3*MU/(3*MU + h))
        beta_1= 1-(S2*S3)
        beta_2= (1-S2)*S3
        deviatoric_tangent_stiffness= (2*MU*beta_1*P_sym) - (2*MU*beta_2*Norm_vector @ np.transpose(Norm_vector))
    k=((3*lamda + 2*MU)/3)
    stress_element=deviatoric_stress + k*trace_strain*Identity_matrixs
    C_tangential=deviatoric_tangent_stiffness + k*(Identity_matrixs@np.transpose(Identity_matrixs))
    A=np.array([2,3,4])
    stress_element=np.delete(stress_element,A,axis=0)
    C_tangential=np.delete(C_tangential,A,axis=0)
    C_tangential=np.delete(C_tangential,A,axis=1)
    for i in range(3,5,1):
        plastic_strain[i][0]=0 # plain strain 
    B_ps[j]=np.transpose(plastic_strain)
    C_al[j][0]=np.transpose(alpha)
    return [C_tangential,stress_element,B_ps,C_al]

def Element_routine(Xe,Element_stiffness_matrixs,MU,lamda,u_element,F_int_Element,Le,thickness_plate,h,B_ps,C_al,P,yield_stress): # Elment routine
    #element routien will return the elements stiffness matrixs
    j=0
    x_values=np.array([-0.57735,-0.57735,0.57735,-0.57735,0.57735,0.57735,-0.57735,0.57735])
    Element_stiffness_matrixs=np.zeros((8,8))
    F_int_Element=np.zeros((8,1))
    #F_int_Element_1=np.zeros((8,1))
    for i in range(0,8,2):
        x1=x_values[i]
        x2=x_values[i+1]
        derivative_x= np.array([[-(1-x2),(1-x2),-(1+x2),(1+x2)],[-(1-x1),-(1+x1),(1-x1),(1+x1)]]) * (1/4)
        Jacobin_matrixs = derivative_x @ Xe
        B_vector= np.linalg.inv(Jacobin_matrixs) @ derivative_x
        B=np.array([[B_vector[0][0],0,B_vector[0][1],0,B_vector[0][2],0,B_vector[0][3],0],[0,B_vector[1][0],0,B_vector[1][1],0,B_vector[1][2],0,B_vector[1][3]],[B_vector[1][0],B_vector[0][0],B_vector[1][1],B_vector[0][1],B_vector[1][2],B_vector[0][2],B_vector[1][3],B_vector[0][3]]])
        [C_tangential,stress_element,B_ps,C_al]=Material_routine(MU,lamda,u_element,B,h,B_ps,C_al,P,j,yield_stress)
        j=j+1 
        Element_stiffness_matrixs= Element_stiffness_matrixs + (np.transpose(B) @ C_tangential @ B)* np.linalg.det(Jacobin_matrixs)*thickness_plate
        F_int_Element= F_int_Element + (np.transpose(B) @ stress_element ) * np.linalg.det(Jacobin_matrixs)*thickness_plate
    #F_int_Element= Element_stiffness_matrixs @ u_element
    #if(np.linalg.norm(F_int_Element)==np.linalg.norm(F_int_Element_1)):
        #print("yes")
    return [Element_stiffness_matrixs,F_int_Element,B_ps,C_al]

#Elastic program for 2D Bilinear Element 
#Here, we are considering a point load as for building the basic strucutre for the required element
#internal parameters
yield_stress=95E6
Youngs_modulus=69E9 #N/meters^2
Poissons_ratio=0.32
h=500e6 #Hardening parameter
alpha=0 #strain like parameter
MU=(Youngs_modulus/(2*(1+Poissons_ratio)))
lamda=((Poissons_ratio*Youngs_modulus)/((1-2*Poissons_ratio)*(1+Poissons_ratio)))
#print(MU)
#print(lamda)
L= eval(input('Enter the length of the plat in meters\n'))
height_plate = eval(input('Enter the height of the plat in meters\n'))
thickness_plate = eval(input("Enter the thickness of the plate meters\n"))
N=eval(input('Enter the number of elements in the x-direction\n')) # No of columns
M=eval(input('Enter the number of elements in the y-direction\n')) # No of rows
Le=L/N #element length
He=height_plate/M #element height
total_nodes= (N+1)*(M+1)
Number_of_elements=N*M
#print(total_nodes)
delta_u=np.zeros((2*total_nodes,1))
Element_stiffness_matrixs=np.zeros((8,8))
Global_F_external=np.zeros((2*total_nodes,1))
Global_plastic_strain=np.zeros((3,M*N))
Global_displacement=np.zeros((2*total_nodes,1))
Global_plastic_strain=np.zeros((Number_of_elements,6))
Gauss_point_plastic_strain=np.zeros((Number_of_elements,4,6))
Gauss_point_alpha=np.zeros((Number_of_elements,4,1))
plastic_strain=np.zeros((6,1))
u_element=np.zeros((8,1))
F_int_Element=np.zeros((8,1))
print("Enter the forces value in newton for each node of interest\n")
print("Enter zero if no forces is applied on the node\n")
# The force part to be made into incremental wise after figuring the flow of the program
for i in range(0,len(Global_F_external),1):
    Global_F_external[i]=eval(input(" "))
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]])
x_disp=np.array([[Le,0],[Le,0],[Le,0],[Le,0]])
y_disp=np.array([[0,He],[0,He],[0,He],[0,He]])
Node_Numbering= np.zeros(((M+1),(N+1)))
#temporary variables for calculations
s=0
s1=0
s2=0
s3=0
s4=0
k=1
k1=0
data={}
for i in range(0,2*total_nodes,2): # creating key values for assembling purpose using even series
    k1=k1+1
    data.update({i:k1})
for i in range(0,M+1): # for segregating the element nodes as a group using a square pattern 
    for j in range(0,N+1):
        Node_Numbering[i][j]=k 
        k=k+1
print(Node_Numbering)

for time_step in range(1,11,1): # time step is split into 10 parts
    R_delta_u=np.ones((2*total_nodes,1))
    Global_displacement=np.zeros((2*total_nodes,1))
    Global_F_external_reduced=(Global_F_external*time_step)/10
    count=0
    while(np.linalg.norm(R_delta_u,np.inf) > (0.005*np.linalg.norm(Global_displacement,np.inf))):
        count=count+1
        Global_stiffness_matrixs=np.zeros((2*total_nodes,2*total_nodes))
        Global_F_internal=np.zeros((2*total_nodes,1))
        P=0
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
                    P=P+1
                B=Gauss_point_plastic_strain[P]
                B=B.reshape((4,6))
                C=Gauss_point_alpha[P]
                C=C.reshape((4,1))
                [Element_stiffness_matrixs,F_int_Element,B,C] = Element_routine(Xe,Element_stiffness_matrixs,MU,lamda,u_element,F_int_Element,Le,thickness_plate,h,B,C,P,yield_stress)
                Global_stiffness_matrixs=Global_stiffness_matrixs + (np.transpose(Assignment_matrixs) @ Element_stiffness_matrixs @ Assignment_matrixs)
                F=(np.transpose(Assignment_matrixs) @ F_int_Element)
                Global_F_internal=Global_F_internal + F
                Gauss_point_plastic_strain[P]=B
                Gauss_point_alpha[P]=C
            if((M*N)>1):
                Xe=Xe + y_disp       
        #print('The obtained stiffness matrixa is:\n')
        #print(Global_stiffness_matrixs)
        #print("Internal_force:\n")
        #print(Global_F_internal)
        G=Global_F_internal-Global_F_external_reduced
        #print("G:",G)
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
        #print(A)
        if(count==1):
            R_delta_u=np.delete(delta_u,A,axis=0)
        Reduced_Global_stiffness_matrix= np.delete(Reduced_Global_stiffness_matrix,A,axis=0)
        Reduced_Global_stiffness_matrix= np.delete(Reduced_Global_stiffness_matrix,A,axis=1)
        Reduced_G=np.delete(Reduced_G,A,axis=0)
        Reduced_displacement=np.delete(Reduced_displacement,A,axis=0)
    #Reduced_F_ext=np.delete(Global_F_external,A,axis=0)
    #calculation part
        #print("reduced k matrixs",Reduced_Global_stiffness_matrix)
        #K_inv=np.linalg.inv(Reduced_Global_stiffness_matrix)
    #R_delta_u= K_inv @ Reduced_F_ext
    #print(R_delta_u)
        #print("Reduced G",Reduced_G)
        R_delta_u=(np.linalg.inv(Reduced_Global_stiffness_matrix)) @ Reduced_G
        #print("r_DELTA_U",R_delta_u)
        Reduced_displacement=Reduced_displacement - R_delta_u
        for i in range(0,len(A),1):
            Reduced_displacement=np.insert(Reduced_displacement,A[i],Global_displacement[A[i]])
        #R_delta_u=np.insert(R_delta_u,A[i],delta_u[i])
        Global_displacement=(Reduced_displacement.reshape(2*total_nodes,1))
    print("displacement:\n",Global_displacement)
    print("Iteration number:",count)