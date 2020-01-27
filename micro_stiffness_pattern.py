import numpy as np 

Element=[2,13,14,16]
print("The flow of dictionary is to be understood")
k=0
nodal_pairing={} # for the local elements indexs num,bering in the k_e system
nodal_i_d={} # to separate the independent and dependent from the local matrixs by assigning 1 or -1.
for i in range(0,8,2):
	nodal_pairing.update({i:Element[k]})
	k=k+1
print(nodal_pairing)

k=1
index_ii={}
for i in range(0,24,2): #### the indexs for ii
	index_ii.update({i:k})
	k=k+1
print(index_ii)
sub_ii=[]
for i in range(0,8,2): # for k_ii index
	for j in range(0,2*len(index_ii),2):
		if(nodal_pairing[i]==index_ii[j]):
			sub_ii.append(j)
			sub_ii.append(j+1)
###For the k_id or k_di
p=13
index_id={}
sub_dd=[]
for i in range(0,8,2):
	index_id.update({i:p})
	p=p+1

for i in range(0,8,2):
	for j in range(0,2*len(index_id),2):
		if(nodal_pairing[i]==index_id[j]):
			sub_dd.append(j)
			sub_dd.append(j+1)
print(index_id)
print("sub_dd",sub_dd)
print(sub_ii)
for i in range(0,8,2):
	if(nodal_pairing[i]<=12):
		nodal_i_d.update({i:1}) ####local element numbering
		nodal_i_d.update({i+1:1})
	else:
		nodal_i_d.update({i:-1})
		nodal_i_d.update({i+1:-1})
#print(nodal_i_d)
p_1=0
p_2=0
p_3=0
p_4=0
q_1=0
q_2=0
q_3=0
q_4=0
k_ii=np.zeros((27,27))
k_id=np.zeros((27,8))
k_di=np.zeros((8,27))
k_dd=np.zeros((8,8))
k=np.ones((8,8))
for i in range(0,8,1): #### finding which part of the k it should go to
	if(q_1==len(sub_ii)):
		p_1=p_1+1
	if(q_2==len(sub_dd)):	
		p_2=p_2+1
	if(q_3==len(sub_ii)):	
		p_3=p_3+1
	if(q_4==len(sub_dd)):	
		p_4=p_4+1
	q_1=0
	q_2=0
	q_3=0
	q_4=0
	for j in range(0,8,1):
		#print(i,j)
		if(nodal_i_d[i]==1) and (nodal_i_d[j]==1):
			#print("The system goes to ind_ind matrix")
			k_ii[sub_ii[p_1]][sub_ii[q_1]]=k[i][j]
			#print(p_1,q_1)
			if(q_1<len(sub_ii)):
				q_1=q_1+1
		elif(nodal_i_d[i]==1) and (nodal_i_d[j]==-1):
			#print("The system goes to ind_dep matrix")
			k_id[sub_ii[p_2]][sub_dd[q_2]]=k[i][j]
			if(q_2<len(sub_dd)):
				q_2=q_2+1
		elif(nodal_i_d[i]==-1)and(nodal_i_d[j]==1):
			#print("The system goes to dep_ind matrix")
			k_di[sub_dd[p_3]][sub_ii[q_3]]=k[i][j]
			if(q_3<len(sub_ii)):
				q_3=q_3+1
		else:
			#print("The system goes to dep_dep matrix")
			k_dd[sub_dd[p_4]][sub_dd[q_4]]=k[i][j]
			if(q_4<len(sub_dd)):
				q_4=q_4+1

print("----------------------------------------------------")
print("k_ii\n")
print(k_ii)
print("----------------------------------------------------")
print("k_id\n")
print(k_id)
print("----------------------------------------------------")
print("k_di\n")
print(k_di)
print("----------------------------------------------------")
print("k_dd\n")
print(k_dd)
print("----------------------------------------------------")
