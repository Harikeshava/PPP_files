import numpy as np
import math as mi
import matplotlib.pyplot as plt
### Micro-scale ###
### The basic theory of the FEM analysis of the micro-scale is that, as I have assumed periodic boundary condition- I have paired the
### boundary nodes. This leads to the concept of independent and dependent nodes wrt strain. The system as a whole has macro strain as 
### the external parameter.For the known displacement of the dependent node, we fidn the displacement of the independent nodes wrt strain.

### The RVE Element is considered with a plate with a hole in the centre
### The micro scale is called at each gauss point of the macro element assumed 

### The below function is for arranging the displacement matrixs for the calculation purpose
def Micro_Displacement_arrangement(u_element,sub_ii,sub_dd,u_independent,u_dependent,nodal_i_d):
	p_1=0
	p_2=0
	for i in range(0,len(u_element),1):
		count_1=0
		count_2=0
		if(nodal_i_d[i]==1):
			count_1=count_1+1
			u_element[i]=u_independent[sub_ii[p_1]]
		else:
			count_2=count_2+1
			u_element[i]=u_dependent[sub_dd[p_2]]
		if(p_1<len(sub_ii)) and (count_1>0):
			p_1=p_1+1
		if(p_2<len(sub_dd)) and (count_2>0):
			p_2=p_2+1
	return [u_element]
### This function is for the calculation of the material model variables.
def Micro_Material_routine(MU,lamda,stress_element,B,u_element,h,B_ps,C_al,P,j,yield_stress):
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
    stress_element=np.delete(stress_element,A,axis=0)
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
    return [C_tangential,stress_element,B_ps,C_al]
###The function for the calculation of the element's variable.
def Micro_Element_routine(x_element,Element_stiffness_micro,F_int_Element,stress_element,u_element,MU,lamda,h,B_ps,C_al,P,yield_stress,thickness_plate):
	j=0
	x_values=np.array([-0.57735,-0.57735,0.57735,-0.57735,-0.57735,0.57735,0.57735,0.57735])
	Element_stiffness_micro=np.zeros((8,8))
	F_int_Element=np.zeros((8,1))
	for k in range(0,8,2):
		x1=x_values[k]
		x2=x_values[k+1]
		derivative_x= np.array([[-(1-x2),(1-x2),-(1+x2),(1+x2)],[-(1-x1),-(1+x1),(1-x1),(1+x1)]]) * (1/4)
		Jacobin_matrixs = derivative_x @ x_element
		B_vector= np.linalg.inv(Jacobin_matrixs) @ derivative_x
		B=np.array([[B_vector[0][0],0,B_vector[0][1],0,B_vector[0][2],0,B_vector[0][3],0],[0,B_vector[1][0],0,B_vector[1][1],0,B_vector[1][2],0,B_vector[1][3]],[B_vector[1][0],B_vector[0][0],B_vector[1][1],B_vector[0][1],B_vector[1][2],B_vector[0][2],B_vector[1][3],B_vector[0][3]]])
		[C_tangential,stress_element,B_ps,C_al]=Micro_Material_routine(MU,lamda,stress_element,B,u_element,h,B_ps,C_al,P,j,yield_stress) 
		j=j+1
		Element_stiffness_micro= Element_stiffness_micro + (np.transpose(B) @ C_tangential @ B)* np.linalg.det(Jacobin_matrixs)*thickness_plate
		F_int_Element= F_int_Element + (np.transpose(B) @ stress_element ) * np.linalg.det(Jacobin_matrixs)*thickness_plate
	return [Element_stiffness_micro,F_int_Element,B_ps,C_al]

def Macro_Material_routine(MU,lamda,B,u_element,Macro_plastic_strain,Macro_drag_alpha,h,yield_stress,j,thickness_plate): # Material routine 
    ############################### Input parameter ########################################
    ### Strain input
    print("The macro strain is given as input now")
    Strain_macro=B @ u_element
    ### units as of now in Meters
    ### Geomentrical variables 
    micro_length=10*10**-6
    micro_height=10*10**-6
    micro_thickness_plate=thickness_plate
    void_radius=2*10**-6
    ### centre of the void in the RVE assumed -(a,b)
    a=(micro_length/2)
    b=(micro_height/2)
    ### Global variables
    global_k_ii=np.zeros((27,27)) #### for the  calculation of the summed k_ii 
    global_k_id=np.zeros((27,8)) #### for the  calculation of the summed k_id
    global_k_di=np.zeros((8,27)) #### for the  calculation of the summed k_di
    global_k_dd=k_dd=np.zeros((8,8)) #### for the  calculation of the summed k_dd
    global_F_int_ii=np.zeros((27,1)) #### for the calculation of the summed f_int_d.
    global_F_int_dd=np.zeros((8,1)) #### for the calculation of the summed f_int_i.
    gauss_point_plastic_strain=Macro_plastic_strain[j]
    gauss_point_alpha=Macro_drag_alpha[j]
    k_total=np.zeros((27,27))#The total stiffness matrixs.
    u_independent=np.zeros((27,1))# The independent element displacement values.
    u_dependent=np.zeros((8,1))# The dependent element displacement values.
    F_int_Element=np.zeros((8,1)) #### for element calculation.
    ### Element variables(local)
    Element_stiffness_micro=np.zeros((8,8))
    F_int_ii=np.zeros((27,1)) ### for the global position arrangement ii.
    F_int_dd=np.zeros((8,1)) ### for the global position arrangement of dd.
    F_int_total=np.zeros((27,1))
    k_ii=np.zeros((27,27)) #The independent-independent stiffness matrixs.
    k_id=np.zeros((27,8))  #The independent-dependent stiffness matrixs.
    k_di=np.zeros((8,27))  #The dependent-independent stiffness matrixs.
    k_dd=np.zeros((8,8))   #The dependent-dependent stiffness matrixs.
    stress_element=np.zeros((3,1)) # local element stress.
    u_element=np.zeros((8,1))# The local element displacements
    ###The co-ordinate variable###
    ##### Pre-processing #####
    ### The pre-processing the stage where, the constructed model is divided into a N number of discrete subregions or element. The elements are connected at the
    ### intersection points called as nodes. The subregions in our analysis represents a geomentry of quadrilateral in shape, these subregions are also called as 
    ### meshes. Now, this developed system's data are used for the  FEM analysis, as an input parameter.
    Micro_coordinate=np.array([[0,0],[(micro_length/2),0],[(void_radius*mi.cos(1.5*mi.pi))+a,(void_radius*mi.sin(1.5*mi.pi))+b],[(void_radius*mi.cos(1.75*mi.pi))+a,(void_radius*mi.sin(1.75*mi.pi))+b],[(void_radius*mi.cos(0*mi.pi))+a,(void_radius*mi.sin(0*mi.pi))+b],[(void_radius*mi.cos(0.25*mi.pi))+a,(void_radius*mi.sin(0.25*mi.pi))+b],[(void_radius*mi.cos(0.5*mi.pi))+a,(void_radius*mi.sin(0.5*mi.pi))+b],[(void_radius*mi.cos(0.75*mi.pi))+a,(void_radius*mi.sin(0.75*mi.pi))+b],[(void_radius*mi.cos(1*mi.pi))+a,(void_radius*mi.sin(1*mi.pi))+b],[(void_radius*mi.cos(1.25*mi.pi))+a,(void_radius*mi.sin(1.25*mi.pi))+b],[0,(micro_height/2)],[0,(micro_height)],[(micro_length/2),micro_height],[micro_length,micro_height],[micro_length,(micro_height/2)],[micro_length,0]])
    print(Micro_coordinate)
    ######## Analysis ########
    ### The data-set obtained from the part of pre-processing is used as the input for our FEM analysis. The code below will solve the system of equation to find the 
    ### unknowns,that is the external parameters. The system has 3 main variables the Total stiffness matrixs, total displacement matrix(the unknowns) and the known 
    ### forces values of each node(external parameter).
    delta_value=np.array([])
    boundary_nodes=[2,1,11,12,13,14,15,16]
    for i in range(0,4,1): ###For caluclating the difference b/w the periodic nodes.
        delta_x=Micro_coordinate[boundary_nodes[i+4]-1][0]-Micro_coordinate[boundary_nodes[i]-1][0]
        delta_y=Micro_coordinate[boundary_nodes[i+4]-1][1]-Micro_coordinate[boundary_nodes[i]-1][1]
        delta_value=np.append(delta_value,delta_x)
        delta_value=np.append(delta_value,delta_y)
    delta_value=delta_value.reshape(8,1)
    ### The below parameters are the global nodal positions of the boundary nodes
    independent_bc=[2,3,0,1,20,21,22,23]
    dependent_bc=[0,1,2,3,4,5,6,7]
    ### This is the matrixs which is equivalent to the assignment matrixs in the macro scale with slight modification
    ### in terms of the periodic boundary conditions
    A_matrix=np.zeros((8,27)) 
    ### The variable k is a temporary scalar variable, used for incremental purpose for the dictionary mapping
    k=0
    ####### Formulation of A_matrixs.
    for i in range(0,8,1):
        A_matrix[dependent_bc[i]][independent_bc[i]]=1
        if(i%2==0) and (i>0):
            k=k+2
        if(i%2==0):
            A_matrix[dependent_bc[i]][24]=delta_value[k]
            A_matrix[dependent_bc[i]][26]=(delta_value[k+1]/2)
        else:
            A_matrix[dependent_bc[i]][25]=delta_value[k+1]
            A_matrix[dependent_bc[i]][26]=(delta_value[k]/2)
    print(A_matrix)
    Micro_element=np.array([[1,2,10,3],[2,16,3,4],[1,10,11,9],[4,16,5,15],[11,9,12,8],[5,15,6,14],[8,7,12,13],[7,6,13,14]])
    k_1=1
    index_ii={} ####The element holds all the index values of independent nodes,in global position
    for i in range(0,24,2): #### the index value of all independent nodes in global matrixs
        index_ii.update({i:k_1})
        k_1=k_1+1
    print(index_ii)
    k_3=13
    index_id={} ###The index of dependent nodes are stored here, the global index values of K_id or K_di
    for i in range(0,8,2):
        index_id.update({i:k_3})
        k_3=k_3+1
    ##### The start of the micro_FEM
    ### The macro strain values are entered here... The coupling part takes place here
    NR_delta_u=np.ones((24,1))
    NR_U=np.zeros((24,1))
    for i in range(0,3,1):
        u_independent[i+24][0]=Strain_macro[i][0]
    print("The element stiffness calculation is going to be started")
    count=0
    ######### Fucntion description ##########
    ###The functions used are:
    ### 1. displacement_arrangement_fn: This function is used for assembling the displacement from global to local variables and vice versa.The dictionary mapping technique
    ### is used here for doing the above function.
    ### 2. Material_routine : This function is used to calculate the solution-dependent state variables of the constitutive mechanical models used for the analysis.
    ### The stress, strains and internal state variables of the model assumed are computed here. The materials current state is kept as the key idea for computing the
    ### the models parameters. Here we use a predictor-corrector analysis for predicting the state of the material, that is weather it is in the elstic state or plastic
    ### state. The von-mises yield criteria is used for the plasticity model.
    ### 3.Element_routine : This function is used to compute the stiffness matrixs and the internal forces acting on the element. This function is called for each element
    ### which is coming in for the analysis. The material routine is called in this function for the calculation purpose, for obtaining the material tangent stiffness 
    ### matrixs and elements stress update(for the calculation of the F_e_interanl).
    ### The element parameters calculation starts here.
    while (np.linalg.norm(NR_delta_u,np.inf) > (0.005*np.linalg.norm(NR_U,np.inf))):
        P=0
        count=count + 1
        global_k_ii=np.zeros((27,27)) 
        global_k_id=np.zeros((27,8))
        global_k_di=np.zeros((8,27))
        global_k_dd=k_dd=np.zeros((8,8))
        global_F_int_ii=np.zeros((27,1))
        global_F_int_dd=np.zeros((8,1))
        #print(np.linalg.norm(NR_delta_u,np.inf))
        #print((0.005*np.linalg.norm(NR_U,np.inf)))
        u_dependent= A_matrix @ u_independent ### The dependent node has relation with strain.
        for i_1 in range(0,8,1):
            print("Element",i_1)
            x_element=np.array([1,1])
            element=Micro_element[i_1] # Elements nodes are assigned to the variable element.
            nodal_pairing={} # for the local elements indexs numbering in the k_e system.
            k_2=0
            for i_2 in range(0,8,2): # dictionary mapping technique
                nodal_pairing.update({i_2:element[k_2]})
                k_2=k_2+1
            print(nodal_pairing)
            sub_ii=[] ### The sub_ii is holding the global positions of the independ nodes coming in for each element, as a local varible for
            #             arrangement purpose in different stiffness matrixs part.
            for i_3 in range(0,8,2): # for k_ii index or k_id or k_di
                for j_1 in range(0,2*len(index_ii),2):
                    if(nodal_pairing[i_3]==index_ii[j_1]):
                        sub_ii.append(j_1)
                        sub_ii.append(j_1+1)
            print(sub_ii)
            sub_dd=[] ### The elemets dependent nodes index are stored here for arragement purpose in K_id or K_di or k_dd
            for i_4 in range(0,8,2): ### Dictionary mapping
                for j_2 in range(0,2*len(index_id),2):
                    if(nodal_pairing[i_4]==index_id[j_2]):
                        sub_dd.append(j_2)
                        sub_dd.append(j_2+1)
            print(sub_dd)
            nodal_i_d={} # to separate the independent and dependent from the local matrixs by assigning 1 or -1,for local k_e
            for i_5 in range(0,8,2):# The loop is used to find which row is of independent and which row is of dependent in local k_e
                if(nodal_pairing[i_5]<=12):
                    nodal_i_d.update({i_5:1}) 
                    nodal_i_d.update({i_5+1:1})
                else:
                    nodal_i_d.update({i_5:-1})
                    nodal_i_d.update({i_5+1:-1})
            print(nodal_i_d)

            for j_3 in range(0,4,1):
                x_element=np.vstack((x_element,Micro_coordinate[element[j_3]-1]))
            x_element=np.delete(x_element,0,axis=0)
            u_element=Micro_Displacement_arrangement(u_element,sub_ii,sub_dd,u_independent,u_dependent,nodal_i_d)
            u_element=np.reshape(u_element,(8,1))
            # temparary local variables
            B=gauss_point_plastic_strain[P]
            B=B.reshape((4,6))
            C=gauss_point_alpha[P]
            C=C.reshape((4,1))
            [Element_stiffness_micro,F_int_Element,B,C]=Micro_Element_routine(x_element,Element_stiffness_micro,F_int_Element,stress_element,u_element,MU,lamda,h,B,C,P,yield_stress,micro_thickness_plate)
            gauss_point_plastic_strain[P]=B
            gauss_point_alpha[P]=C
            P=P+1
            ### Arrangment of the dependent and independent nodes stiffness values in 4 different combination matixs.
            ### The below scalar variables,before for loop,are the temporary variables for calculation.
            ### The scalar variables
            p_1=0
            p_2=0
            p_3=0
            p_4=0
            q_1=0
            q_2=0
            q_3=0
            q_4=0
            ### The matrix variables
            k_ii=np.zeros((27,27)) #The independent-independent stiffness matrixs.
            k_id=np.zeros((27,8))  #The independent-dependent stiffness matrixs.
            k_di=np.zeros((8,27))  #The dependent-independent stiffness matrixs.
            k_dd=np.zeros((8,8))   #The dependent-dependent stiffness matrixs.
            for i_6 in range(0,8,1): #### finding which part of the k it should go to
                if(q_1==len(sub_ii)) and (len(sub_ii)>0):
                    p_1=p_1+1
                if(q_2==len(sub_dd)) and (len(sub_dd)>0):	
                    p_2=p_2+1
                if(q_3==len(sub_ii)) and (len(sub_ii)>0):	
                    p_3=p_3+1
                if(q_4==len(sub_dd)) and (len(sub_dd)>0):	
                    p_4=p_4+1
                q_1=0
                q_2=0
                q_3=0
                q_4=0
                for j_4 in range(0,8,1):
                    #print(i,j)
                    if(nodal_i_d[i_6]==1) and (nodal_i_d[j_4]==1):
                        #print("The system goes to ind_ind matrix")
                        k_ii[sub_ii[p_1]][sub_ii[q_1]]= Element_stiffness_micro[i_6][j_4]
                        #print(p_1,q_1)
                        if(q_1<len(sub_ii)):
                            q_1=q_1+1
                    elif(nodal_i_d[i_6]==1) and (nodal_i_d[j_4]==-1):
                        #print("The system goes to ind_dep matrix")
                        k_id[sub_ii[p_2]][sub_dd[q_2]]=Element_stiffness_micro[i_6][j_4]
                        if(q_2<len(sub_dd)):
                            q_2=q_2+1
                    elif(nodal_i_d[i_6]==-1)and(nodal_i_d[j_4]==1):
                        #print("The system goes to dep_ind matrix")
                        k_di[sub_dd[p_3]][sub_ii[q_3]]=Element_stiffness_micro[i_6][j_4]
                        if(q_3<len(sub_ii)):
                            q_3=q_3+1
                    else:
                        #print("The system goes to dep_dep matrix")
                        k_dd[sub_dd[p_4]][sub_dd[q_4]]=Element_stiffness_micro[i_6][j_4]
                        if(q_4<len(sub_dd)):
                            q_4=q_4+1
            p_1=0
            p_2=0
            F_int_ii=np.zeros((27,1)) 
            F_int_dd=np.zeros((8,1))
            for i_7 in range(0,8,1):
                count_1=0
                count_2=0
                if(nodal_i_d[i_7]==1):
                    count_1=count_1+1
                    F_int_ii[sub_ii[p_1]][0]=F_int_Element[i_7]
                else:
                    count_2=count_2+1
                    F_int_dd[sub_dd[p_2]][0]=F_int_Element[i_7]
                if(p_1<len(sub_ii))and(count_1>0):
                    p_1 = p_1+1
                if(p_2<len(sub_dd))and(count_2>0):
                    p_2 = p_2+1
            global_F_int_ii= global_F_int_ii + F_int_ii
            global_F_int_dd= global_F_int_dd + F_int_dd
            print("----------------------------------------------------")
            print("k_ii\n")
            #print(k_ii)
            global_k_ii= global_k_ii +k_ii
            print("----------------------------------------------------")
            print("k_id\n")
            #print(k_id)
            global_k_id=global_k_id + k_id
            print("----------------------------------------------------")
            print("k_di\n")
            #print(k_di)
            global_k_di=global_k_di + k_di
            print("----------------------------------------------------")
            print("k_dd\n")
            #print(k_dd)
            global_k_dd= global_k_dd + k_dd
            print("----------------------------------------------------")
            print("#########")
            print("Calculation of k_total\n")
        #k_total_1=np.zeros((27,27))
        #k_total_1= k_ii + (k_id @ A_matrix) + (np.transpose(A_matrix) @ k_di ) + (np.transpose(A_matrix) @ k_dd @ A_matrix)
        k_total= global_k_ii + (global_k_id @ A_matrix) + (np.transpose(A_matrix) @ global_k_di ) + (np.transpose(A_matrix) @ global_k_dd @ A_matrix)
        #print(np.allclose(k_total,k_total_1))
        F= np.transpose(A_matrix) @ global_F_int_dd
        F_int_total= global_F_int_ii + np.transpose(A_matrix) @ global_F_int_dd
        print("K_total:\n")
        #print(k_total)
        print("----------------------------------------------------")
        print("F_int_total:")
        #print(F_int_total)
        print("----------------------------------------------------")
        print("The newton-raphsons method is going to start here\n")
        NR_stiffness_matrix= k_total[:24,:24]
        A=[16,17]
        print("The stiffness matrixs of 24*24 is formed\n")
        print("The stationary node 9's nodal values are removed")
        NR_stiffness_matrix=np.delete(NR_stiffness_matrix,A,axis=0)
        NR_stiffness_matrix=np.delete(NR_stiffness_matrix,A,axis=1)
        NR_G= F_int_total[:24]
        NR_G=np.delete(NR_G,A,axis=0)
        NR_U= u_independent[:24]
        NR_U= np.delete(NR_U,A,axis=0)
        NR_delta_u=(np.linalg.inv(NR_stiffness_matrix)) @ NR_G
        NR_U= NR_U - NR_delta_u
        for i in range(0,len(A),1):
            NR_U=np.insert(NR_U,A[i],0)
        NR_U=(NR_U.reshape(24,1))
        for i in range(0,len(NR_U),1):
            u_independent[i][0]=NR_U[i][0]
    print("count:\n",count)
    print("u_independent\n",u_independent)
    print("u_dependent\n",u_dependent)
    volume= micro_length * micro_height * micro_thickness_plate
    Macro_stress_element= (F[24:,0])/volume
    Macro_stress_element=(Macro_stress_element.reshape(3,1))
    k_3_3=k_total[24:,24:] 
    k_3_24= k_total[24:,:24]
    k_inv= k_total[:24,:24]
    k_24_3= k_total[:24,24:]
    term=k_3_24 @ np.linalg.inv(k_inv) @ k_24_3
    Macro_C_tangential= (k_3_3 - term)/volume 
    print(Macro_C_tangential)
    return [Macro_C_tangential,Macro_stress_element,Macro_plastic_strain,Macro_drag_alpha]

##### Macro scale program #####
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

def Macro_Element_routine(Xe,Element_stiffness_matrixs,MU,lamda,u_element,F_int_Element,thickness_plate,h,Macro_plastic_strain,Macro_drag_alpha,yield_stress): # Element routine
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
        [C_tangential,stress_element,Macro_plastic_strain,Macro_drag_alpha]=Macro_Material_routine(MU,lamda,B,u_element,Macro_plastic_strain,Macro_drag_alpha,h,yield_stress,j,thickness_plate)
        j=j+1 
        Element_stiffness_matrixs= Element_stiffness_matrixs + (np.transpose(B) @ C_tangential @ B)* np.linalg.det(Jacobin_matrixs)*thickness_plate
        F_int_Element= F_int_Element + (np.transpose(B) @ stress_element ) * np.linalg.det(Jacobin_matrixs)*thickness_plate
    #F_int_Element= Element_stiffness_matrixs @ u_element
    #if(np.linalg.norm(F_int_Element)==np.linalg.norm(F_int_Element_1)):
        #print("yes")
    return [Element_stiffness_matrixs,F_int_Element,Macro_plastic_strain,Macro_drag_alpha]

#Elastic program for 2D Bilinear Element- as the homogenisation technique is applied
#Here, we are considering a point load as for building the basic strucutre for the required element
############################ Macro Input parameter #################################
### units as of now in Meters ###
### The variables that can me changed for analysis ###
### Geomentrical variables
macro_length=100 * 10**-2 #(100cm)
macro_height=50 * 10**-2  #(50cm)
thickness_plate=500*10**-2  #(500cm)
### Material property variables
Youngs_modulus= 70E9 #N/meter^2
Poissons_ratio=0.30
yield_stress=95E6 #N/meter^2
h=200e6 #Hardening parameter
force_mag = 100 # N/m^2 # force distribution  magnitude
N=1 #eval(input('number of elements in the x-direction\n')) # No of columns
M=1 #eval(input('number of elements in the y-direction\n')) # No of rows
#####################################################################################


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
L= macro_length 
height_plate = macro_height 
thickness_plate = thickness_plate
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
Gauss_point_plastic_strain=np.zeros((Number_of_elements,4,8,4,6))
Gauss_point_alpha=np.zeros((Number_of_elements,4,8,4,1))
u_element=np.zeros((8,1))
F_int_Element=np.zeros((8,1))
#print("Enter the forces value in newton for each node of interest\n")
#print("Enter zero if no forces is applied on the node\n")
## Meshing Variables ##
Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]])
x_disp=np.array([[Le,0],[Le,0],[Le,0],[Le,0]])
y_disp=np.array([[0,He],[0,He],[0,He],[0,He]])
##############################################
############ Testing Parameters ##############
# For verifying the test cases, follow the instructions in the manual section (to run the program for the test cases).
#tension_disp=np.array([[0,0],[-2.110*10**-7,-6.588*10**-7],[0,0],[1.85*10**-7,-6.28*10**-7]])
#Xe=np.array([[0,0],[Le,0],[0,He],[Le,He]]) + tension_disp
#Global_displacement[2][0]= 4.909*10**-7 
#Global_displacement[3][0]= 7.748*10**-8
#Global_displacement[4][0]= -4.909*10**-7
#Global_displacement[5][0]= -7.748*10**-8
#### Boundary condition ####
#A=np.array([[2],[3],[4],[5]]) The variable A array should be uncommented in the line 638 of the program line below, where the array A holds the boundary condiitons.
##############################################
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
## Force Initialization ## For runnin the program for test cases the below 2 lines should be commented 
#Global_F_external[2][0]= force_mag
#Global_F_external[6][0]= force_mag
for j in range(0,M+1,1):    
    for n, a in data.items():    # for printing the value's key
        if a == Node_Numbering[j][N]:
            i=n
            Global_F_external[i][0]= force_mag
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

for time_step in range(5,105,5): # 5 % increment at each iteration
    R_delta_u=np.ones((2*total_nodes,1))
    Global_displacement=np.zeros((2*total_nodes,1))
    Global_F_external_reduced=(0.5 * macro_height*thickness_plate*Global_F_external*time_step)/100
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
                Macro_plastic_strain=Gauss_point_plastic_strain[P]
                Macro_drag_alpha=Gauss_point_alpha[P]
                [Element_stiffness_matrixs,F_int_Element,Macro_plastic_strain,Macro_drag_alpha] = Macro_Element_routine(Xe,Element_stiffness_matrixs,MU,lamda,u_element,F_int_Element,thickness_plate,h,Macro_plastic_strain,Macro_drag_alpha,yield_stress)
                Global_stiffness_matrixs=Global_stiffness_matrixs + (np.transpose(Assignment_matrixs) @ Element_stiffness_matrixs @ Assignment_matrixs)
                F=(np.transpose(Assignment_matrixs) @ F_int_Element)
                Global_F_internal=Global_F_internal + F
                Gauss_point_plastic_strain[P]=Macro_plastic_strain
                Gauss_point_alpha[P]=Macro_drag_alpha
            if((M*N)>1):
                Xe=Xe + y_disp       
        #print('The obtained stiffness matrixa is:\n')
        #print(Global_stiffness_matrixs)
        #print("Internal_force:\n")
        #print(Global_F_internal)
        G=Global_F_internal-Global_F_external_reduced
        #print("G:",G)
        Reduced_Global_stiffness_matrix=Global_stiffness_matrixs
        X=Global_stiffness_matrixs -np.transpose(Global_stiffness_matrixs)
        print(X)
        print(np.linalg.det(Global_stiffness_matrixs))
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
        #### Boundary condition ####
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
plt.plot((Final_position[0])*10**10,(Final_position[1])*10**10)
plt.xlabel("x(m)")
plt.ylabel("y(m)")
plt.title("Final positions of the node")
plt.show()

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

#The stress and strain plots along x and y
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


fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, strain_x)
lx = plt.xlabel("x (m)")
ly = plt.ylabel("y (m)")
plt.title("strain xx")
fig.colorbar(cmap)
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