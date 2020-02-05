import numpy as np 
import math as mi 
import matplotlib.pyplot as plt
#### The below code is developed for the analysis of the FEM analysis for a plate with a hole in Macro-scale using the Quadrilateral element.

### The function for arrangement of the u_element values from the global displacement values
def displacement_arrangement_fn(Global_node_position,Global_displacement,s,u_element,t):
    for n, a in Global_node_position.items():    # for printing the value's key
        if a == s:
            i=n
            u_element[t]=Global_displacement[i]
            u_element[t+1]=Global_displacement[i+1]
    return u_element
### The function for assignment matrixs arrangement.
def assignment_matrix_fn(Global_node_position,s,A_matrix,t):
    for n, a in Global_node_position.items():    # for printing the value's key
        if a == s:
            i=n
            A_matrix[t][i]=1
            A_matrix[t+1][i+1]=1
    return A_matrix
### The function for the calculation of the material model variables.
def Material_routine(MU,lamda,B,u_element,H,B_ps,C_al,P,j,yield_stress,stress_33):
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
    drag_stress=H*alpha
    plastic_fn=mag_deviatoric-(2/3)**(1/2)*(yield_stress + drag_stress)
    if(plastic_fn<10**-6):
        print("Elastic calculation")
        deviatoric_stress=deviatoric_stress_trial
        deviatoric_tangent_stiffness= 2*MU*P_sym
        plastic_strain=np.transpose(plastic_strain)
        alpha=np.transpose(alpha)
    else:
        print("Plastic calculation")
        plastic_corrector= (plastic_fn)/(2*MU + (2/3)*H)
        deviatoric_stress = deviatoric_stress_trial - 2*MU*Norm_vector*plastic_corrector
        plastic_strain= np.transpose(plastic_strain) + plastic_corrector*Norm_vector
        alpha = np.transpose(alpha) + (2/3)**(1/2)*(plastic_corrector)
        S2 = (plastic_fn/mag_deviatoric)
        S3 = (3*MU/(3*MU + H))
        beta_1= 1-(S2*S3)
        beta_2= (1-S2)*S3
        deviatoric_tangent_stiffness= (2*MU*beta_1*P_sym) - (2*MU*beta_2*Norm_vector @ np.transpose(Norm_vector))
    k=((3*lamda + 2*MU)/3)
    stress_element=deviatoric_stress + k*trace_strain*Identity_matrixs
    C_tangential=deviatoric_tangent_stiffness + k*(Identity_matrixs@np.transpose(Identity_matrixs))
    A=np.array([2,3,4])
    stress_33[j][0]=stress_element[2][0]
    stress_element=np.delete(stress_element,A,axis=0)
    C_tangential=np.delete(C_tangential,A,axis=0)
    C_tangential=np.delete(C_tangential,A,axis=1)
    for i in range(3,5,1):
        plastic_strain[i][0]=0 # plain strain 
    B_ps[j]=np.transpose(plastic_strain)
    C_al[j][0]=np.transpose(alpha)
    return [C_tangential,stress_element]
###The function for the calculation of the element's variable.
def Element_routine(x_element,Element_stiffness_matrix_macro,F_internal_Element,u_element,MU,lamda,H,B_ps,C_al,P,yield_stress,stress_33,thickness_plate):
    j=0
    x_values=np.array([-0.57735,-0.57735,0.57735,-0.57735,0.57735,0.57735,-0.57735,0.57735])
    Element_stiffness_matrix_macro=np.zeros((8,8))
    F_internal_Element=np.zeros((8,1))
    for k in range(0,8,2):
        x1=x_values[k]
        x2=x_values[k+1]
        derivative_x= np.array([[-(1-x2),(1-x2),-(1+x2),(1+x2)],[-(1-x1),-(1+x1),(1-x1),(1+x1)]]) * (1/4)
        Jacobin_matrixs = derivative_x @ x_element
        B_vector= np.linalg.inv(Jacobin_matrixs) @ derivative_x
        B=np.array([[B_vector[0][0],0,B_vector[0][1],0,B_vector[0][2],0,B_vector[0][3],0],[0,B_vector[1][0],0,B_vector[1][1],0,B_vector[1][2],0,B_vector[1][3]],[B_vector[1][0],B_vector[0][0],B_vector[1][1],B_vector[0][1],B_vector[1][2],B_vector[0][2],B_vector[1][3],B_vector[0][3]]])
        [C_tangential,stress_element]=Material_routine(MU,lamda,B,u_element,H,B_ps,C_al,P,j,yield_stress,stress_33) 
        Element_stiffness_matrix_macro= Element_stiffness_matrix_macro + (np.transpose(B) @ C_tangential @ B)* np.linalg.det(Jacobin_matrixs)*thickness_plate
        j=j+1
        F_internal_Element= F_internal_Element + (np.transpose(B) @ stress_element ) * np.linalg.det(Jacobin_matrixs)*thickness_plate
    return [Element_stiffness_matrix_macro,F_internal_Element]


############################ Input parameter #################################
### units as of now in Meters ###
### Geomentrical variables
macro_length=100 * 10**-2 #(100cm)
macro_height=50 * 10**-2  #(50cm)
thickness_plate=1*10**-2  #(1cm)
void_radius=2*10**-2      #(2cm)
n=8 ### Number of point in the circumfrenece of the circle.
### centre of the circle-(a,b)
a=(macro_length/2)
b=(macro_height/2)
### Material property variables
Youngs_modulus=210E9 #N/meters
Poissons_ratio=0.30
yield_stress=95E6
H=500e6 #Hardening parameter
MU=(Youngs_modulus/(2*(1+Poissons_ratio)))
lamda=((Poissons_ratio*Youngs_modulus)/((1-2*Poissons_ratio)*(1+Poissons_ratio)))
### Elements and node variables
total_nodes=16
Number_of_elements=8
Global_node_position={}
#Macro_element=np.array([[1,2,10,3],[2,16,3,4],[1,10,11,9],[4,16,5,15],[11,9,12,8],[5,15,6,14],[8,7,12,13],[7,6,13,14]])
Macro_element=np.array([[1,2,10,11],[2,3,11,12],[3,4,12,13],[4,5,13,14],[5,6,14,15],[6,7,15,16],[7,8,16,9],[8,1,9,10]]) ### The element numbering
### Global variables
Global_displacement=np.zeros((2*total_nodes,1))
Global_stiffness_Matrix=np.zeros((2*total_nodes,2*total_nodes))
Global_F_internal=np.zeros((2*total_nodes,1))
Global_F_external=np.zeros((2*total_nodes,1))
Global_displacement=np.zeros((2*total_nodes,1))
Global_plastic_strain=np.zeros((Number_of_elements,6)) #### got to add this for my convergences of the system later...
Gauss_point_plastic_strain=np.zeros((Number_of_elements,4,6))
Gauss_point_alpha=np.zeros((Number_of_elements,4,1))
Gauss_point_stress_33=np.zeros((Number_of_elements,4,1)) #### got to check this stress_33 storage part.
### Local  variables for the incoming element for analysis
u_element=np.zeros((8,1))
Element_stiffness_matrix_macro=np.zeros((8,8))
stress_33=np.zeros((4,1))
plastic_strain=np.zeros((6,1))
F_internal_Element=np.zeros((8,1)) #### Force _element_internal
local_positions=[0,2,4,6] ### The local posiitons
### Newton-Raphsons parameters
delta_u=np.zeros((2*total_nodes,1))
###Reduced Newton-Raphson parameters
R_delta_u=np.ones((2*total_nodes,1))


##### Pre-processing #####
### The pre-processing the stage where, the constructed model is divided into a N number of discrete subregions or element. The elements are connected at the
### intersection points called as nodes. The subregions in our analysis represents a geomentry of quadrilateral in shape, these subregions are also called as 
### meshes. Now, this developed system's data are used for the  FEM analysis, as an input parameter.
print("The macro-coordinates are calculated here\n")
### Meshing of the given system into quadrilateral element is done here
Macro_coordinate=np.array([[0,0],[(macro_length/2),0],[macro_length,0],[macro_length,(macro_height/2)],[macro_length,macro_height],[(macro_length/2),macro_height],[0,macro_height],[0,(macro_height/2)]])
for i in range(4,12,1):
    x= mi.cos(2*np.pi/n*i)*void_radius + a
    y= mi.sin(2*np.pi/n*i)*void_radius + b
    c=np.array([x,y])
    c=np.reshape(c,(1,2))
    Macro_coordinate=np.append(Macro_coordinate,c,axis=0)
print(Macro_coordinate)
# The varibale k is used as a temporary varibal in the below for loop, for setting up  the dictionary mapping techniques in the othe parts of the code.
k=1
for i in range(0,(2*total_nodes),2):
    Global_node_position.update({i:k})
    k=k+1
print("Global_node_positions\n",Global_node_position)
print("The local_node _positions\n",local_positions)

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
print("The External forces are entered now\n")
print("Enter zero if no forces is applied on the node\n")
for i in range(0,len(Global_F_external),1): # Point Load of interest
   Global_F_external[i]=eval(input(" "))
for time_step in range(1,11,1): ### The iteration is split into 10 steps
    R_delta_u=np.ones((2*total_nodes,1))
    Global_displacement=np.zeros((2*total_nodes,1))
    Global_F_external_reduced=(Global_F_external*time_step)/10
    count=0
    while (np.linalg.norm(R_delta_u,np.inf) > (0.005*np.linalg.norm(Global_displacement,np.inf))):
        count=count+1
        Global_stiffness_Matrix=np.zeros((2*total_nodes,2*total_nodes))
        Global_F_internal=np.zeros((2*total_nodes,1))
        P=0
        for i in range(0,Number_of_elements,1): ### loop for element
            A_matrix=np.zeros((8,2*total_nodes))  ### the total nodes is fixed as 16
            print("Element",i+1)
            x_element=np.array([1,1])
            element=Macro_element[i] # Elements nodes are assigned to the variable element.
            for j in range(0,4,1):
                x_element=np.vstack((x_element,Macro_coordinate[element[j]-1]))
            x_element=np.delete(x_element,0,axis=0)
            print("Element nodes:\n",element)
            print("The global co-ordinate values\n",x_element)
            ### The elements displacement assignment
            dummy=np.arange(32)
            for i in range(0,4,1):
                u_element=displacement_arrangement_fn(Global_node_position,Global_displacement,element[i],u_element,local_positions[i])
            ### The elements assignment matrixs arrangment
            for i in range(0,4,1):
                A_matrix=assignment_matrix_fn(Global_node_position,element[i],A_matrix,local_positions[i])
            ### The variables B,C,D,P are temporary variables for calculation purpose
            B=Gauss_point_plastic_strain[P]
            B=B.reshape((4,6))
            C=Gauss_point_alpha[P]
            C=C.reshape((4,1))
            D=(Gauss_point_stress_33[P])
            stress_33=D.reshape((4,1))
            [Element_stiffness_matrix_macro,F_internal_Element]=Element_routine(x_element,Element_stiffness_matrix_macro,F_internal_Element,u_element,MU,lamda,H,B,C,P,yield_stress,stress_33,thickness_plate)
            Global_stiffness_Matrix=Global_stiffness_Matrix + (np.transpose(A_matrix) @ Element_stiffness_matrix_macro @ A_matrix)
            Global_F_internal=Global_F_internal + (np.transpose(A_matrix) @ F_internal_Element)
            Gauss_point_plastic_strain[P]=B
            Gauss_point_alpha[P]=C
            Gauss_point_stress_33[P]=stress_33
            P=P+1
        print("Global-stiffness-micro:\n",Global_stiffness_Matrix)
        X=Global_stiffness_Matrix - np.transpose(Global_stiffness_Matrix)
        print(X)
        print("Global-F-internal:\n",Global_F_internal)
        #### Newton-Raphsons method
        G_vector=Global_F_internal-Global_F_external_reduced 
        Reduced_Global_stiffness_matrix=Global_stiffness_Matrix
        Reduced_displacement=Global_displacement
        Reduced_G=G_vector
        A=[] 
        #Reduction of matirxs sizes
        Node_Numbering=np.array([[1],[7],[8]])
        for j in range(0,3,1):    # got to change this 
            for n, a in Global_node_position.items():    # for printing the value's key
                if a == Node_Numbering[j][0]:
                    i=n
                    A.append(i)
                    A.append(i+1)
        A=np.asarray(A)
        if(count==1):
            R_delta_u=np.delete(delta_u,A,axis=0)
        Reduced_Global_stiffness_matrix= np.delete(Reduced_Global_stiffness_matrix,A,axis=0)
        Reduced_Global_stiffness_matrix= np.delete(Reduced_Global_stiffness_matrix,A,axis=1)
        Reduced_G=np.delete(Reduced_G,A,axis=0)
        Reduced_displacement=np.delete(Reduced_displacement,A,axis=0)
        R_delta_u=(np.linalg.inv(Reduced_Global_stiffness_matrix)) @ Reduced_G
        Reduced_displacement=Reduced_displacement - R_delta_u
        for i in range(0,len(A),1):
            Reduced_displacement=np.insert(Reduced_displacement,A[i],Global_displacement[A[i]])
        Global_displacement=(Reduced_displacement.reshape(2*total_nodes,1))
    print("displacement:\n",Global_displacement)
    print("Iteration number:",count)    
#plt.plot(Micro_coordinate[:,0],Micro_coordinate[:,1])
#plt.show()