import numpy as np 
import math as mi 
import matplotlib.pyplot as plt
### The code is for the micro-scale FEM analysis
### The basic theory of the FEM analysis of the micro-scale is that, as we have assumed periodic boundary condition- me have paired the
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
def Micro_Material_routine(MU,lamda,stress_element,B,u_element,h,B_ps,C_al,P,j,yield_stress,stress_33):
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
	stress_33[j][0]= stress_element[2][0] #lamda*(Strain_element[0][0]+Strain_element[1][0]) - 2*MU*plastic_strain[2][0]
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
def Micro_Element_routine(x_element,Element_stiffness_micro,F_int_Element,stress_element,u_element,MU,lamda,h,B_ps,C_al,P,yield_stress,stress_33,thickness_plate):
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
		[C_tangential,stress_element,B_ps,C_al]=Micro_Material_routine(MU,lamda,stress_element,B,u_element,h,B_ps,C_al,P,j,yield_stress,stress_33) 
		j=j+1
		Element_stiffness_micro= Element_stiffness_micro + (np.transpose(B) @ C_tangential @ B)* np.linalg.det(Jacobin_matrixs)*thickness_plate
		F_int_Element= F_int_Element + (np.transpose(B) @ stress_element ) * np.linalg.det(Jacobin_matrixs)*thickness_plate
	return [Element_stiffness_micro,F_int_Element,B_ps,C_al]

############################### Input parameter ########################################
### units as of now in Meters
### Geomentrical variables
micro_length=10*10**-6
micro_height=10*10**-6
thickness_plate=500*10**-2
void_radius=2*10**-6
### Material Properties ###
Youngs_modulus=70E9 #N/meters
Poissons_ratio=0.30
yield_stress=95E6
h=200e6 #Hardening parameter
Strain_macro=np.array([[5.16754*10**-7],[-1.0335*10**-7],[-4.268*10**-7]])
###########################################################################################
MU=(Youngs_modulus/(2*(1+Poissons_ratio)))
lamda=((Poissons_ratio*Youngs_modulus)/((1-2*Poissons_ratio)*(1+Poissons_ratio)))
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
Gauss_point_plastic_strain=np.zeros((8,4,6))
Gauss_point_alpha=np.zeros((8,4,1))
Gauss_point_stress_33=np.zeros((8,4,1)) 
### Element variables(local)
Element_stiffness_micro=np.zeros((8,8))
F_int_ii=np.zeros((27,1)) ### for the global position arrangement ii.
F_int_dd=np.zeros((8,1)) ### for the global position arrangement of dd.
F_int_Element=np.zeros((8,1)) #### for element calculation.
F_int_total=np.zeros((27,1))
k_ii=np.zeros((27,27)) #The independent-independent stiffness matrixs.
k_id=np.zeros((27,8))  #The independent-dependent stiffness matrixs.
k_di=np.zeros((8,27))  #The dependent-independent stiffness matrixs.
k_dd=np.zeros((8,8))   #The dependent-dependent stiffness matrixs.
k_total=np.zeros((27,27))#The total stiffness matrixs.
u_independent=np.zeros((27,1))# The independent element displacement values.
u_dependent=np.zeros((8,1))# The dependent element displacement values.
stress_element=np.zeros((3,1)) # local element stress.
u_element=np.zeros((8,1))# The local element displacements
###The co-ordinate variable###
##### Pre-processing #####
### The pre-processing the stage where, the constructed model is divided into a N number of discrete subregions or element. The elements are connected at the
### intersection points called as nodes. The subregions in our analysis represents a geomentry of quadrilateral in shape, these subregions are also called as 
### meshes. Now, this developed system's data are used for the  FEM analysis, as an input parameter.
Micro_coordinate=np.array([[0,0],[(micro_length/2),0],[(void_radius*mi.cos(1.5*mi.pi))+a,(void_radius*mi.sin(1.5*mi.pi))+b],[(void_radius*mi.cos(1.75*mi.pi))+a,(void_radius*mi.sin(1.75*mi.pi))+b],[(void_radius*mi.cos(0*mi.pi))+a,(void_radius*mi.sin(0*mi.pi))+b],[(void_radius*mi.cos(0.25*mi.pi))+a,(void_radius*mi.sin(0.25*mi.pi))+b],[(void_radius*mi.cos(0.5*mi.pi))+a,(void_radius*mi.sin(0.5*mi.pi))+b],[(void_radius*mi.cos(0.75*mi.pi))+a,(void_radius*mi.sin(0.75*mi.pi))+b],[(void_radius*mi.cos(1*mi.pi))+a,(void_radius*mi.sin(1*mi.pi))+b],[(void_radius*mi.cos(1.25*mi.pi))+a,(void_radius*mi.sin(1.25*mi.pi))+b],[0,(micro_height/2)],[0,(micro_height)],[(micro_length/2),micro_height],[micro_length,micro_height],[micro_length,(micro_height/2)],[micro_length,0]])
print(Micro_coordinate)
plt.plot(Micro_coordinate[:,0],Micro_coordinate[:,1])
plt.show()
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
### Strain input
print("The macro strain is given as input now")
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
	print(np.linalg.norm(NR_delta_u,np.inf))
	print((0.005*np.linalg.norm(NR_U,np.inf)))
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
		B=Gauss_point_plastic_strain[P]
		B=B.reshape((4,6))
		C=Gauss_point_alpha[P]
		C=C.reshape((4,1))
		D=(Gauss_point_stress_33[P])
		stress_33=D.reshape((4,1))
		[Element_stiffness_micro,F_int_Element,B,C]=Micro_Element_routine(x_element,Element_stiffness_micro,F_int_Element,stress_element,u_element,MU,lamda,h,B,C,P,yield_stress,stress_33,thickness_plate)
		Gauss_point_plastic_strain[P]=B
		Gauss_point_alpha[P]=C
		#Gauss_point_stress_33[P]=stress_33
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
		print(k_ii)
		global_k_ii= global_k_ii +k_ii
		print("----------------------------------------------------")
		print("k_id\n")
		print(k_id)
		global_k_id=global_k_id + k_id
		print("----------------------------------------------------")
		print("k_di\n")
		print(k_di)
		global_k_di=global_k_di + k_di
		print("----------------------------------------------------")
		print("k_dd\n")
		print(k_dd)
		global_k_dd= global_k_dd + k_dd
		print("----------------------------------------------------")
		print("#########")
		print("Calculation of k_total\n")
	#k_total_1=np.zeros((27,27))
	#k_total_1= k_ii + (k_id @ A_matrix) + (np.transpose(A_matrix) @ k_di ) + (np.transpose(A_matrix) @ k_dd @ A_matrix)
	k_total= global_k_ii + (global_k_id @ A_matrix) + (np.transpose(A_matrix) @ global_k_di ) + (np.transpose(A_matrix) @ global_k_dd @ A_matrix)
	#print(np.allclose(k_total,k_total_1))
	F_int_total= global_F_int_ii + np.transpose(A_matrix) @ global_F_int_dd
	print("K_total:\n")
	print(k_total)
	print("----------------------------------------------------")
	print("F_int_total:")
	print(F_int_total)
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
volume= micro_length * micro_height * thickness_plate
#stress_element= (F[24:,0])/volume
k_3_3=k_total[24:,24:] 
k_3_24= k_total[24:,:24]
k_inv= k_total[:24,:24]
k_24_3= k_total[:24,24:]
term=k_3_24 @ np.linalg.inv(k_inv) @ k_24_3
t_1= (k_3_3 /volume )
C_tangential= t_1 - ((term)/volume)
print(C_tangential)