import numpy as np 
import itertools as it
import matplotlib.pyplot as plt
from  matplotlib import cm
def u_arrangement(data,Global_displacement,s,u_element,t): # The function for arranging  the displacement element accessed
    for n, a in data.items():    # for printing the value's key
        if a == s:
            i=n
            u_element[t][0]=Global_displacement[i][0]
            u_element[t+1][0]=Global_displacement[i+1][0]
    return u_element

def assemble_function(data,Assignment_matrixs,s,j): # The function to assemble the assignment matrixs
    for n, a in data.items():    # for printing the value's key
        if a == s:
            i=n
            Assignment_matrixs[j][i]=1
            Assignment_matrixs[j+1][i+1]=1
    return(Assignment_matrixs)

def Material_routine(MU,lamda,u_element,B,h,B_ps,C_al,P,j,yield_stress,stress_33): # Material routine 
    #material routine will return the C_t matrixs,stress at gauss points and the internal state variables
    Norm_vector=np.zeros((6,1))
    P_sym=np.array([[(2/3),(-1/3),(-1/3),0,0,0],[(-1/3),(2/3),(-1/3),0,0,0],[(-1/3),(-1/3),(2/3),0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    Strain_element=np.zeros((3,1))
    Strain_element= B @ u_element
    for i in range(2,5,1):
        Strain_element= np.insert(Strain_element,i,0,axis=0)
    for i in range(3,6,1):
       Strain_element[i][0]=(Strain_element[i][0])/2
    print("Strain_element:\n",Strain_element)
    trace_strain= Strain_element[0][0]+Strain_element[1][0]+Strain_element[2][0]
    Identity_matrixs=np.array([[1],[1],[1],[0],[0],[0]])
    deviatoric_strain= Strain_element - (trace_strain*Identity_matrixs)/3
    plastic_strain=(B_ps[j].reshape(1,6))
    #print(plastic_strain)
    alpha=C_al[j][0]
    deviatoric_stress = 2*MU*(deviatoric_strain-np.transpose(plastic_strain))
    #deviatoric_stress[2][0]=-2*MU*plastic_strain[0][2] ### replacing the strain 33 component
    deviatoric_stress_trial=deviatoric_stress
    s=0
    s1=1
    for i in range(0,6,1):
        if(i>2):
            s1=2
        s=s+(deviatoric_stress_trial[i][0]*deviatoric_stress_trial[i][0])*s1
    mag_deviatoric=(s)**(1/2)
    print("Von Mises:",mag_deviatoric)
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
        plastic_corrector= ((plastic_fn)/(2*MU + (2*h)/3))
        deviatoric_stress = deviatoric_stress_trial - 2*MU*Norm_vector*plastic_corrector
        s=0
        s1=1
        for i in range(0,6,1):
            if(i>2):
                s1=2
            s=s+(deviatoric_stress[i][0]*deviatoric_stress[i][0])*s1
        mag_deviatoric_1=(s)**(1/2)
        plastic_strain= np.transpose(plastic_strain) + plastic_corrector*Norm_vector
        #plastic_strain[2][0]=-(plastic_strain[0][0]+plastic_strain[1][0])
        alpha = np.transpose(alpha) + (2/3)**(1/2)*(plastic_corrector)
        S2 = (plastic_fn/mag_deviatoric)
        S3 = (3*MU/(3*MU + h))
        beta_1= 1-(S2*S3)
        beta_2= (1-S2)*S3
        term_1=Norm_vector @ np.transpose(Norm_vector)
        deviatoric_tangent_stiffness= (2*MU*beta_1*P_sym) - (2*MU*beta_2*term_1)
    k= lamda + (2/3)*MU 
    stress_element=deviatoric_stress + k*trace_strain*Identity_matrixs
    term = k * (Identity_matrixs @ np.transpose(Identity_matrixs))
    C_tangential=deviatoric_tangent_stiffness + term
    A=np.array([2,3,4])
    stress_33[j][0]= stress_element[2][0] #lamda*(Strain_element[0][0]+Strain_element[1][0]) - 2*MU*plastic_strain[2][0]
    stress_element=np.delete(stress_element,A,axis=0)
    print("Stress_element",stress_element)
    C_tangential=np.delete(C_tangential,A,axis=0)
    C_tangential=np.delete(C_tangential,A,axis=1)
    C_tangential[0][2]=0
    C_tangential[1][2]=0
    C_tangential[2][0]=0
    C_tangential[2][1]=0
    C_tangential[2][2]= (C_tangential[2][2]/2)
    for i in range(3,5,1):
        plastic_strain[i][0]=0 # plain strain 
    B_ps[j]=np.transpose(plastic_strain)
    C_al[j][0]=np.transpose(alpha)
    return [C_tangential,stress_element,B_ps,C_al,stress_33]

def Element_routine(Xe,Element_stiffness_matrixs,MU,lamda,u_element,F_int_Element,Le,thickness_plate,h,B_ps,C_al,P,yield_stress,stress_33): # Element routine
    #element routien will return the elements stiffness matrixs
    j=0
    x_values=np.array([-0.57735,-0.57735,0.57735,-0.57735,-0.57735,0.57735,0.57735,0.57735])
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
        [C_tangential,stress_element,B_ps,C_al,stress_33]=Material_routine(MU,lamda,u_element,B,h,B_ps,C_al,P,j,yield_stress,stress_33)
        j=j+1 
        Element_stiffness_matrixs= Element_stiffness_matrixs + (np.transpose(B) @ C_tangential @ B)* np.linalg.det(Jacobin_matrixs)*thickness_plate
        F_int_Element= F_int_Element + (np.transpose(B) @ stress_element ) * np.linalg.det(Jacobin_matrixs)*thickness_plate
    #F_int_Element= Element_stiffness_matrixs @ u_element
    #if(np.linalg.norm(F_int_Element)==np.linalg.norm(F_int_Element_1)):
        #print("yes")
    return [Element_stiffness_matrixs,F_int_Element,B_ps,C_al,stress_33]

#Elastic program for 2D Bilinear Element- as the homogenisation technique is applied
#Here, we are considering a point load as for building the basic strucutre for the required element
############################ Input parameter #################################
### units as of now in Meters ###
### The variables that can me changed for analysis ###
### Geomentrical variables
macro_length=100 * 10**-2 #(100cm)
macro_height=50 * 10**-2  #(50cm)
thickness_plate=500*10**-2  #(500cm)
### Material property variables
Youngs_modulus= 70E9 #N/meters
Poissons_ratio=0.30
yield_stress=95E6
h=200e6 #Hardening parameter
force_mag = 300*10**5 # N/m^2 # forces distributions magnitude
N=1 #eval(input('number of elements in the x-direction\n')) # No of columns
M=1 #eval(input('number of elements in the y-direction\n')) # No of rows
###############################################################################
MU=(Youngs_modulus/(2*(1+Poissons_ratio)))
lamda=((Poissons_ratio*Youngs_modulus)/((1-2*Poissons_ratio)*(1+Poissons_ratio)))
### Local  variables for the incoming element for analysis
u_element=np.zeros((8,1))
stress_33=np.zeros((4,1))
plastic_strain=np.zeros((6,1))
F_internal_Element=np.zeros((8,1)) #### Force _element_internal
##### Pre-processing #####
### The pre-processing the stage where, the constructed model is divided into a N number of discrete subregions or element. The elements are connected at the
### intersection points called as nodes. The subregions in our analysis represents a geomentry of quadrilateral in shape, these subregions are also called as 
### meshes. Now, this developed system's data are used for the  FEM analysis, as an input parameter.
L= macro_length #eval(input('Enter the length of the plat in meters\n'))
height_plate = macro_height #eval(input('Enter the height of the plat in meters\n'))
thickness_plate = thickness_plate #eval(input("Enter the thickness of the plate meters\n"))
### Elements and node variables
Le=L/N #element length
He=height_plate/M #element height
total_nodes= (N+1)*(M+1) # the total number of nodes in the system
Number_of_elements=N*M
#print(total_nodes)
### Global variables
delta_u=np.zeros((2*total_nodes,1))
Element_stiffness_matrixs=np.zeros((8,8))
Global_F_external=np.zeros((2*total_nodes,1))
Global_displacement=np.zeros((2*total_nodes,1))
Gauss_point_plastic_strain=np.zeros((Number_of_elements,4,6))
Gauss_point_alpha=np.zeros((Number_of_elements,4,1))
Gauss_point_stress_33=np.zeros((Number_of_elements,4,1)) 
stress_33=np.zeros((4,1))
plastic_strain=np.zeros((6,1))
u_element=np.zeros((8,1))
F_int_Element=np.zeros((8,1))
## The forces initialization ## ### For running the program for test cases the below 2 lines should be commented 
Global_F_external[2][0]= force_mag
Global_F_external[6][0]= force_mag
#### Meshing variables ####
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]])
x_disp=np.array([[Le,0],[Le,0],[Le,0],[Le,0]])
y_disp=np.array([[0,He],[0,He],[0,He],[0,He]])
############################
######## Testing Parameters ######
#For verifying the test cases, follow the instructions in the manual section (to run the program for the test cases)
#tension_disp=np.array([[0,0],[5.55*10**-9,1.88*10**-10],[0,0],[-6.33*10**-9,1.97*10**-10]])
#Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) #+ tension_disp
#Global_displacement[2][0]= 4.909*10**-7 
#Global_displacement[3][0]= 7.748*10**-8
#Global_displacement[4][0]= -4.909*10**-7
#Global_displacement[5][0]= -7.748*10**-8
#### Boundary condition ####
#A=np.array([[2],[3],[4],[5]]) The variable A array should be uncommented in the line 293 of the program line below, where the array A holds the boundary condiitons
#############################
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
for i in range(0,2*total_nodes,2): # The variable data consists of the global node positions, the even series of the numbers are arranged here.
    k1=k1+1
    data.update({i:k1})
for i in range(0,M+1): # for segregating the element nodes as a group using a square pattern 
    for j in range(0,N+1):
        Node_Numbering[i][j]=k 
        k=k+1
print(Node_Numbering)
######## Analysis ########
### The data-set obtained from the part of pre-processing is used as the input for our FEM analysis. The code below will solve the system of equation to find the 
### unknowns,that is the external parameters. The system has 3 main variables the Total stiffness matrixs, total displacement matrix(the unknowns) and the known 
### forces values of each node(external parameter).
######### Fucntion description ##########
###The functions used are:
### 1. displacement_arrangement_fn: This function is used for assembling the displacement from global to local variables and vice versa.The dictionary mapping technique
### is used here for doing the above function.
### 2. assignment_matrix_fn : This function is used for the assembling  the assignment matrixs for the element coming in as the local element.The assigment matrixs
### is used for the mapping the element stiffness matrix and global stiffness matrix.
### 3. Material_routine : This function is used to calculate the solution-dependent state variables of the constitutive mechanical models used for the analysis.
### The stress, strains and internal state variables of the model assumed are computed here. The materials current state is kept as the key idea for computing the
### the models parameters. Here we use a predictor-corrector analysis for predicting the state of the material, that is weather it is in the elstic state or plastic
### state. The von-mises yield criteria is used for the plasticity model.
### 4.Element_routine : This function is used to compute the stiffness matrixs and the internal forces acting on the element. This function is called for each element
### which is coming in for the analysis. The material routine is called in this function for the calculation purpose, for obtaining the material tangent stiffness 
### matrixs and elements stress update(for the calculation of the F_e_interanl).
x=[]
y=[]
for time_step in range(1,11,1): # time step is split into 10 parts
    R_delta_u=np.ones((2*total_nodes,1))
    Global_displacement=np.zeros((2*total_nodes,1)) ## comment this line for making it displacement driven
    Global_F_external_reduced=(thickness_plate*Global_F_external*time_step)/10
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
                D=(Gauss_point_stress_33[P])
                stress_33=D.reshape((4,1))
                [Element_stiffness_matrixs,F_int_Element,B,C,stress_33] = Element_routine(Xe,Element_stiffness_matrixs,MU,lamda,u_element,F_int_Element,Le,thickness_plate,h,B,C,P,yield_stress,stress_33)
                Global_stiffness_matrixs=Global_stiffness_matrixs + (np.transpose(Assignment_matrixs) @ Element_stiffness_matrixs @ Assignment_matrixs)
                F=(np.transpose(Assignment_matrixs) @ F_int_Element)
                Global_F_internal=Global_F_internal + F
                Gauss_point_plastic_strain[P]=B
                Gauss_point_alpha[P]=C
                Gauss_point_stress_33[P]=stress_33
            if((M*N)>1):
                Xe=Xe + y_disp       
        #print('The obtained stiffness matrixa is:\n')
        #print(Global_stiffness_matrixs)
        print("Internal_force:\n")
        print(Global_F_internal)
        G=Global_F_internal-Global_F_external_reduced
        #print("G:",G)
        Reduced_Global_stiffness_matrix=Global_stiffness_matrixs
        X=Global_stiffness_matrixs -np.transpose(Global_stiffness_matrixs)
        print(X)
        Reduced_displacement=Global_displacement
        Reduced_G=G
        A=[] #### Boundary condition ####
        #Reduction of matirxs sizes
        for j in range(0,M+1,1):    
            for n, a in data.items():    # for printing the value's key
                if a == Node_Numbering[j][0]:
                    i=n
                    A.append(i)
                    A.append(i+1)
        A=np.asarray(A) #### Boundary condition ####
        #print(A)
        #### The below A array values should be used for the test cases
        #A=np.array([[2],[3],[4],[5]])
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

x=[]
y=[]
for i in range(0,len(Global_displacement),1):
    if(i%2==0):
        x.append(Global_displacement[i][0])
    else:
        y.append(Global_displacement[i][0])
print(x,y)
Final_displacment = np.append(np.vstack(x),np.vstack(y),axis=1)
Final_position= Xe + Final_displacment
plt.scatter(Xe[:,0]*10**10,Xe[:,1]*10**10)
plt.scatter((Final_position[:,0]*10**10),(Final_position[:,1]*10**10))
plt.xlabel("x(m)")
plt.ylabel("y(m)")
plt.title("Final positions of the node")
plt.show()

print("The virtual work of the system with irregular meshing is:\n")
print( (Global_stiffness_matrixs@ Global_displacement) - Global_F_external)
#### Post-processing ####
Element_x=np.linspace(0,macro_length,50)
Element_y=np.linspace(0,macro_height,50)
U_x=np.zeros((len(Element_x),len(Element_y)))
U_y=np.zeros((len(Element_x),len(Element_y)))
Aera= macro_length * macro_height
xx , yy = np.meshgrid(Element_x,Element_y)
# interpolation of the displacements using the obtained final displacement from the analysis
for i in range(0,len(xx),1):
    for j in range(0,len(yy),1):
        N_1=  ((xx[i][j]-macro_length)*(yy[i][j]-macro_height))*(1/Aera)
        N_2=  - ((xx[i][j])*(yy[i][j]-macro_height))*(1/Aera)
        N_3=  - ((xx[i][j]-macro_length)*(yy[i][j]))*(1/Aera)
        N_4=  ((xx[i][j])*(yy[i][j]))*(1/Aera)
        U_x[i][j]= N_1*x[0] + N_2 * x[1] + N_3 * x[2] + N_4 * x[3]
        U_y[i][j]= N_1 * y[0] + N_2 * y[1] + N_3 * y[2] + N_4 * y[3]
print(U_x)
print(U_y)
# The displacement plots along x and y direction
fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, U_x)
lx = plt.xlabel("x(m)")
ly = plt.ylabel("y(m)")
plt.title("Displacement of the nodes along x-direction")
fig.colorbar(cmap)
plt.show(fig)
fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, U_y)
lx = plt.xlabel("x (m)")
ly = plt.ylabel("y (m)")
plt.title("Displacement of the nodes along y-direction")
fig.colorbar(cmap)
plt.show(fig)
# The stress interpolations
strain_x=np.zeros((len(Element_x),len(Element_y)))
strain_y=np.zeros((len(Element_x),len(Element_y)))
strain_xy=np.zeros((len(Element_x),len(Element_y)))
stress_x=np.zeros((len(Element_x),len(Element_y)))
stress_y=np.zeros((len(Element_x),len(Element_y)))
stress_xy=np.zeros((len(Element_x),len(Element_y)))
stress_vector=np.zeros((3,1))
C_tangential=np.array([[2*MU+lamda,lamda,0],[lamda,2*MU+lamda,0],[0,0,MU]])
for i in range(0,len(xx),1):
    for j in range(0,len(yy),1):
        # Partial derivative w.r.t x
        N_1_x=  ((yy[i][j]-macro_height))*(1/Aera)
        N_2_x=  - ((yy[i][j]-macro_height))*(1/Aera)
        N_3_x=  - ((yy[i][j]))*(1/Aera)
        N_4_x=  ((yy[i][j]))*(1/Aera)
        # Partial derivative w.r.t y
        N_1_y=  ((xx[i][j]-macro_length))*(1/Aera)
        N_2_y=  - ((xx[i][j]))*(1/Aera)
        N_3_y=  - ((xx[i][j]-macro_length))*(1/Aera)
        N_4_y=  ((xx[i][j]))*(1/Aera)
        strain_x[i][j]= N_1_x*x[0] + N_2_x * x[1] + N_3_x * x[2] + N_4_x * x[3]
        strain_y[i][j]= N_1_y * y[0] + N_2_y * y[1] + N_3_y * y[2] + N_4_y * y[3]
        strain_xy[i][j]= N_1_x*y[0] + N_2_x * y[1] + N_3_x * y[2] + N_4_x * y[3] + N_1_y * x[0] + N_2_y * x[1] + N_3_y * x[2] + N_4_y * x[3]
        strain=np.array([[strain_x[i][j]],[strain_y[i][j]],[strain_xy[i][j]]])
        stress_vector = C_tangential @ strain
        stress_x[i][j]= stress_vector[0][0]
        stress_y[i][j]= stress_vector[1][0]
        stress_xy[i][j]= stress_vector[2][0]

print(strain_x)
fig, ax = plt.subplots()
cmap = ax.contourf(xx, yy, strain_x)
lx = plt.xlabel("x (m)")
ly = plt.ylabel("y (m)")
plt.title("strain xx")
plt.colorbar(cmap)
plt.show(fig)

fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, strain_y)
lx = plt.xlabel("x (m)")
ly = plt.ylabel("y (m)")
plt.title("strain yy")
fig.colorbar(cmap)
plt.show(fig)

fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, strain_xy)
lx = plt.xlabel("x (m)")
ly = plt.ylabel("y (m)")
plt.title("strain xy")
fig.colorbar(cmap)
plt.show(fig)

## stress
fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, stress_x)
lx = plt.xlabel("x (m)")
ly = plt.ylabel("y (m)")
plt.title("stress xx (N/m^2)")
fig.colorbar(cmap)
plt.show(fig)

fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, stress_y)
lx = plt.xlabel("x (m)")
ly = plt.ylabel("y (m)")
plt.title("stress yy (N/m^2)")
fig.colorbar(cmap)
plt.show(fig)

fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, stress_xy)
lx = plt.xlabel("x (m)")
ly = plt.ylabel("y (m)")
plt.title("stress xy (N/m^2)")
fig.colorbar(cmap)
plt.show(fig)