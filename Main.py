import math 
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

G = 6.673e-11 #Gravity Const.
M1 = 1.9891e30 #Sun Mass (kg)
M2 = 5.9723e24 # Earth mass (kg)
T = 365.25*24*3660 #s day * hours *seconds
R = 1.4710e11 #Perihellion
Q = (G*M1*T**2)/R**3
#r_stateVector = np.array([1.4710e11, 0]) *1000   # m, z= 4207.43331
#v_StateVector = np.array([0, 30290]) #m/s ,z=3.679500667
#a_acceleration = np.array([-C*r_stateVector[0] / linalg.norm(r_stateVector)**3,
                            #(-C*r_stateVector[1] / linalg.norm(r_stateVector)**3)*(T/R)])


a = 1.4960e11 #semi-major axis
e = .0167  #eccentricity
n = math.sqrt(G*M1*M2/a**3)
p = a*(1-e**2)


Mean = []
def MeanAnomaly(k,h):
    tp = 0                               #time at perigee
    M = 0
    for i in range(0,k):
        M = M+ n*(h-tp)
        Mean.append(M)
    return M
Z = MeanAnomaly(10000, .0001)

#print(Mean)

Enn = [] 
L2 = [] 
   
def NewtonsMethod(E, k, tol):
    def BC(x):
        AA = 1-e*math.cos(E)
        return AA    
    counter = 1
    En = E
    for i in range(0, k):
        def AB(E):
            BB = E - e*math.sin(E) - Mean[i]
            return BB
        if BC(E) == 0:
            return E
        En = E - AB(E)/BC(E)
        if abs(En-E) < tol:
            #print(i, En)
            Enn.append(En)
            L2.append(AB(En))
            break
        else:
            E = En
            Enn.append(E)
            L2.append(AB(En))
            counter += 1
            #print(i, E)

    return En
En = NewtonsMethod( 0, 10000, .0001)

RO = []      
def rPolarCoord(k):
    for i in range(k):        
        rO = a*(1-e*math.cos(Enn[i]))
        RO.append(rO)
    return rO
rPOLAR = rPolarCoord(10000)
TA = []
def TrueAnom(k):
    for i in range(k):
        phi = 180*(math.acos((a*(1-e**2)-RO[i])/(e*RO[i]))) / math.pi
        TA.append(phi)
    return phi
TrueAn = TrueAnom(10000)
#print(RO)
#print(TA)

XX = []
YY = []  
def Cart(k):
    for i in range(0,k):    
        X = a*(math.cos(Enn[i])-e)
        Y = a*(math.sqrt(1-e**2))*math.sin(Enn[i])
        XX.append(X)
        YY.append(Y)
    return X
    return Y
Cartesian = Cart(10000)
#print(XX)
#print(YY)



fig = plt.figure()
ax = plt.axes(projection='3d')
xline = [XX]
yline = [YY]
zline = 0
plt.title('Exact')
ax.scatter3D(xline, yline, zline, color='pink')
plt.show()



X_POS = 1 # as X_POS = xpos / R
Y_POS = 0 
VX_POS = 0
VY_POS = 30290*(T/R)

  
XPRunge = np.zeros(100001)
YPRunge = np.zeros(100001)
VxRunge = np.zeros(100001)
VyRunge = np.zeros(100001)

XPRunge[0] = X_POS

VyRunge[0] = VY_POS

h = .0001
k = 10000


XPOS =np.zeros(10001)
YPOS =np.zeros(10001)
VELX =np.zeros(10001)
VELY =np.zeros(10001)

XPOS[0] = X_POS

VELY[0] = VY_POS

def RungeKutta2(x,vy,k):    
    def XPrime(x, y, vx, vy):
        X = vx
        return X
    def YPrime(x, y, vx, vy):
        Y = vy
        return Y
    def VxPrime(x, y, vx, vy):
        VX = (-Q * x) / (math.sqrt(x**2 + y**2))**3
        return VX
    def VyPrime(x, y, vx, vy):
        VY = (-Q * y) / (math.sqrt(x**2 + y**2))**3
        return VY 
    for i in range(k-1):
        KX1 = h*XPrime(XPOS[i],YPOS[i],VELX[i],VELY[i])
        KY1 = h*YPrime(XPOS[i],YPOS[i],VELX[i],VELY[i])
        KVX1 = h*VxPrime(XPOS[i],YPOS[i],VELX[i],VELY[i])
        KVY1 = h*VyPrime(XPOS[i],YPOS[i],VELX[i],VELY[i])
        KX2 = h* XPrime(XPOS[i]+KX1 /2,YPOS[i]+KY1 /2,VELX[i]+KVX1 /2,VELY[i]+KVY1 /2)
        KY2 =  h* YPrime(XPOS[i]+KX1 /2,YPOS[i]+KY1 /2,VELX[i]+KVX1 /2,VELY[i]+KVY1 /2)
        KVX2 = h* VxPrime(XPOS[i]+KX1 /2,YPOS[i]+KY1 /2,VELX[i]+KVX1 /2,VELY[i]+KVY1 /2)
        KVY2 = h* VyPrime(XPOS[i]+KX1 /2,YPOS[i]+KY1 /2,VELX[i]+KVX1 /2,VELY[i]+KVY1 /2)
        XPOS[i+1] = XPOS[i] + KX2
        YPOS[i+1] = YPOS[i] + KY2 
        VELX[i+1] = VELX[i] + KVX2
        VELY[i+1] = VELY[i] + KVY2
        
    return XPOS
    return YPOS
    return VELX
    return VELY

RK2 = RungeKutta2(X_POS,VY_POS ,10000)   

plt.figure()
ax1 = plt.axes(projection='3d')
xline1 = [XPOS]
yline1 = [YPOS]
zline1 = 0
plt.title('2nd Order Runge Kutta')
ax1.scatter3D(xline1, yline1, zline1, color='red')
plt.show() 
        
def RungeKutta4(x,vy,k): 
    #XPRunge = np.zeros(k)
    #YPRunge = np.zeros(k)
    #VxRunge = np.zeros(k)
    #VyRunge = np.zeros(k)   
    def XPrime(x, y, vx, vy):
        X = vx
        return X
    def YPrime(x, y, vx, vy):
        Y = vy
        return Y
    def VxPrime(x, y, vx, vy):
        VX = (-Q * x) / (math.sqrt(x**2 + y**2))**3
        return VX
    def VyPrime(x, y, vx, vy):
        VY = (-Q * y) / (math.sqrt(x**2 + y**2))**3
        return VY 
    for i in range(k-1):
        Kx1 = h*XPrime(XPRunge[i],YPRunge[i],VxRunge[i],VyRunge[i])
        Ky1 = h*YPrime(XPRunge[i],YPRunge[i],VxRunge[i],VyRunge[i])
        Kvx1 = h*VxPrime(XPRunge[i],YPRunge[i],VxRunge[i],VyRunge[i])
        Kvy1 = h*VyPrime(XPRunge[i],YPRunge[i],VxRunge[i],VyRunge[i])
        RKx2 = h*XPrime(XPRunge[i]+(Kx1)/2,YPRunge[i]+(Ky1)/2,VxRunge[i]+(Kvx1)/2,VyRunge[i]+(Kvy1)/2)
        RKy2 =h*YPrime(XPRunge[i]+(Kx1)/2,YPRunge[i]+(Ky1)/2,VxRunge[i]+(Kvx1)/2,VyRunge[i]+(Kvy1)/2)
        RKvx2 =h*VxPrime(XPRunge[i]+(Kx1)/2,YPRunge[i]+(Ky1)/2,VxRunge[i]+(Kvx1)/2,VyRunge[i]+(Kvy1)/2)
        RKvy2 =h*VyPrime(XPRunge[i]+(Kx1)/2,YPRunge[i]+(Ky1)/2,VxRunge[i]+(Kvx1)/2,VyRunge[i]+(Kvy1)/2)  
        Kx3 = h*XPrime(XPRunge[i]+(RKx2)/2,YPRunge[i]+(RKy2)/2,VxRunge[i]+(RKvx2)/2,VyRunge[i]+(RKvy2)/2)
        Ky3 = h*YPrime(XPRunge[i]+(RKx2)/2,YPRunge[i]+(RKy2)/2,VxRunge[i]+(RKvx2)/2,VyRunge[i]+(RKvy2)/2)
        Kvx3 = h*VxPrime(XPRunge[i]+(RKx2)/2,YPRunge[i]+(RKy2)/2,VxRunge[i]+(RKvx2)/2,VyRunge[i]+(RKvy2)/2)
        Kvy3 = h*VyPrime(XPRunge[i]+(RKx2)/2,YPRunge[i]+(RKy2)/2,VxRunge[i]+(RKvx2)/2,VyRunge[i]+(RKvy2)/2)
        RKx4 = h*XPrime(XPRunge[i]+Kx3,YPRunge[i]+Ky3,VxRunge[i]+Kvx3,VyRunge[i]+Kvy3)
        RKy4 = h*YPrime(XPRunge[i]+Kx3,YPRunge[i]+Ky3,VxRunge[i]+Kvx3,VyRunge[i]+Kvy3)
        RKvx4 = h*VxPrime(XPRunge[i]+Kx3,YPRunge[i]+Ky3,VxRunge[i]+Kvx3,VyRunge[i]+Kvy3)
        RKvy4 = h*VyPrime(XPRunge[i]+Kx3,YPRunge[i]+Ky3,VxRunge[i]+Kvx3,VyRunge[i]+Kvy3)
        XPRunge[i+1] = XPRunge[i] + (1/6)*(Kx1+(2*RKx2)+(2*Kx3)+RKx4)
        YPRunge[i+1] = YPRunge[i] + (1/6)*(Ky1+(2*RKy2)+(2*Ky3)+RKy4)
        VxRunge[i+1] = VxRunge[i] + (1/6)*(Kvx1+(2*RKvx2)+(2*Kvx3)+RKvx4)
        VyRunge[i+1] = VyRunge[i] + (1/6)*(Kvy1+(2*RKvy2)+(2*Kvy3)+RKvy4)

    return XPRunge
    return YPRunge
    return VxRunge
    return VyRunge
    

RK = RungeKutta4(X_POS,VY_POS ,10000)
    

plt.figure()
ax1 = plt.axes(projection='3d')
xline1 = [XPRunge]
yline1 = [YPRunge]
zline1 = 0
plt.title('4th Order Runge Kutta')
ax1.scatter3D(xline1, yline1, zline1, color='green')
plt.show() 


#print(XX)
#print(YY)

Energy = -(G*M1*M2) / (2*a)
#print(Energy)


Energy2 = np.zeros(9999)
def EnergyRK2(k):
    for i in range(k-1):
       X =XPOS[i+1]
       Y =YPOS[i+1] 
       VX = VELX[i+1] 
       VY = VELY[i+1] 
       Energy2[i] = .5 * (VX**2 + VY**2) -(Q /(math.sqrt(X**2 + Y**2)))
    return Energy
LL = EnergyRK2(10000)
#print(Energy2)

Energy1 = np.zeros(9999)
def EnergyRK4(k):
    for i in range(k-1):
        X = XPRunge[i+1] 
        Y = YPRunge[i+1] 
        VX = VxRunge[i+1] 
        VY = VyRunge[i+1] 
        Energy1[i] = .5 * (VX**2 + VY**2) -Q /(math.sqrt(X**2 + Y**2))
    return Energy
L = EnergyRK4(10000)
#print(Energy1)

plt.figure()
plt.plot(Energy1)
plt.title('4th Order Runge Kutta Energy')
plt.show()


plt.figure()
plt.plot(Energy2)
plt.title('2nd Order Runge Kutta Energy')
plt.show()

v = np.zeros(10000)
def VISVIVA(k):
    for i in range(k-1):
        v[i] = math.sqrt(G*M1*M2((2/RO[i]) - (1/a)))
    return v
plt.figure()
plt.plot(v)
plt.title('Exact Energy')
plt.show()

for i in range(10000):
    print('Exact Solution:', i, math.sqrt(XX[i]**2 + YY[i]**2) * 6.68459e-12)
    print('Runge Kutta 2:', i, math.sqrt(XPOS[i]**2 + YPOS[i]**2))
    print('Runge Kutta 4:', i, math.sqrt(XPRunge[i]**2 + YPRunge[i]**2))
