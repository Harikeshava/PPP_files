import numpy as np 
import matplotlib.pyplot as plt
length=10
height=50
x=np.linspace(0,length,50)
y=np.linspace(0,height,50)
U_x=np.zeros((len(x),len(y)))
U_y=np.zeros((len(x),len(y)))
Ones_value=np.ones((len(x),len(y)))
Aera= length * height
x_disp=[0,4.9,0,4.9]
y_disp=[0,0.774,0,-0.774]
xx , yy = np.meshgrid(x,y)
for i in range(0,len(xx),1):
    for j in range(0,len(yy),1):
        N_1=  ((xx[i][j]-length)*(yy[i][j]-height))*(1/Aera)
        N_2=  - ((xx[i][j])*(yy[i][j]-height))*(1/Aera)
        N_3=  - ((xx[i][j]-length)*(yy[i][j]))*(1/Aera)
        N_4=  ((xx[i][j])*(yy[i][j]))*(1/Aera)
        U_x[i][j]= N_1*x_disp[0] + N_2 * x_disp[1] + N_3 * x_disp[2] + N_4 * x_disp[3]
        U_y[i][j]= N_1 * y_disp[0] + N_2 * y_disp[1] + N_3 * y_disp[2] + N_4 * y_disp[3]
print(U_x)
print(U_y)

fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, U_x)
lx = plt.xlabel("x(*10^-8 m)")
ly = plt.ylabel("y(*10^-8 m)")
plt.title("Displacement of the nodes along x-direction")
fig.colorbar(cmap)
plt.show(fig)
fig, ax = plt.subplots()
cmap = plt.contourf(xx, yy, U_y)
lx = plt.xlabel("x (*10^-8 m)")
ly = plt.ylabel("y (*10^-8 m)")
plt.title("Displacement of the nodes along y-direction")
fig.colorbar(cmap)
plt.show(fig)