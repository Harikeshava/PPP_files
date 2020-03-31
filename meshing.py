import numpy as np 
import math 
import matplotlib.pyplot as plt
def PointsInCircum(r,n,a,b):
    X=[]
    Y=[]
    for i in range(0,n+1):
        x= math.cos(2*np.pi/n*i)*r + a
        print(x)
        X.append(x)
        y= math.sin(2*np.pi/n*i)*r + b
        Y.append(y)
    return [X,Y]
print("The mesh generation is done here,for a plate with a hole")
Length=eval(input("Enter the length of the recangle"))
height=eval(input("Enter the height of the rectangle"))
radius_1=eval(input("Enter the radius of the circle/hole"))
n=eval(input("Enter the number of point that you need in the circumfrence of the circle"))
a=Length/2
b=height/2
distance=((Length-0)**2 +(height-0)**2)**(1/2)
radius_2=int(distance/2)
X=[]
Y=[]
s=n-2
[x,y]=PointsInCircum(radius_1,s,a,b)
plt.scatter(x,y)
for i in range(int(radius_1)+1,int(radius_2),1):
    [x,y]=PointsInCircum(i,s,a,b)
    for j in range(0,len(x),1):
       if(((x[j]>=0) and (x[j]<=Length)) and ((y[j]>=0) and (y[j]<=height))):
            X.append(x[j])
            Y.append(y[j])
x_r=np.linspace(0,Length,int(n/4)+1)
y_r=np.linspace(0,height,int(n/4)+1)
x_c_1=[0]*(int(n/4)+1)
x_c_2=[Length]*(int(n/4)+1)
y_c_1=[0]*(int(n/4)+1)
y_c_2=[height]*(int(n/4)+1)
plt.scatter(x_r,y_c_1)
plt.scatter(x_c_2,y_r)
plt.scatter(x_r,y_c_2)
plt.scatter(x_c_1,y_r)
plt.scatter(X,Y)
plt.xlabel("Length of the RVE")
plt.ylabel("Height of the RVE")
plt.show()