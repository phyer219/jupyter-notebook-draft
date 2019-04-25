import matplotlib.pyplot as plt
import math

n=1000 #the number of steps
nt=10
dt=3*math.pi/nt
y1=[0]*(n+1)
y2=[0]*(n+1)
y=[0]*2
y1[0]=y[0]=0 #initial values
y2[0]=y[1]=2

def rungeKutta(y,t,dt):
    c1=[0]*2
    c2=[0]*2
    c3=[0]*2
    c4=[0]*2

    c1=f(y,t)
    for i in range(2):
        c2[i]=y[i]+dt*c1[i]/2
    c2=f(c2,t+dt/2)
    for i in range(2):
        c3[i]=y[i]+dt*c2[i]/2
    c3=f(c3,t+dt/2)
    for i in range(2):
        c4[i]=y[i]+dt*c3[i]
    c4=f(c4,t+dt)
    for i in range(2):
        c1[i]=y[i]+dt*(c1[i]+2*c2[i]+2*c3[i]+c4[i])/6
    return c1

def f(y,t):
    g=0.5 #the two parameters of the model
    b=0.6

    v=[0]*2 #velocity
    v[0]=y[1]
    v[1]=-g*y[1]-y[0]**3+b*math.cos(t)
    return v #return f(1),f(2)

for i in range(n): #main
        t=dt*i
        y=rungeKutta(y,t,dt)
        y1[i+1]=y[0]
        y2[i+1]=y[1]

x1=[]
x2=[]
for i in range(0,n):
    x1.append(y1[i])
    x2.append(y2[i])

plt. scatter (x1,x2)
plt. show ()

    
    

