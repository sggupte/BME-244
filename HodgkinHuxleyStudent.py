import numpy as np
import matplotlib.pyplot as plt

# Define all the initial values for the simulation
gNa=    #gNa max (1/ohm)                         
gK=     #gK max (1/ohm)         
GO=     #Go (1/ohm)               
Cm=     #Cell membrane capacitance
Vr=     #Resting Membrane Potential (V)
ENa=    #Nernst Potential for Sodium (V)
EK=     #Nernst Potential of Potassium (V)
EO=     #Nernst Potential of all other ions (V)
                             
# Define all the gating variables
M0=     #M gates (no units)                     
N0=     #N gates (no units)
H0=     #H gates (no units)  

# Stimulus current applied between 7 and 7.2 (in Amps)
Istim=                    

# Start time, end time, and time increment (in seconds)
t0=0
tEnd=40
dt=0.01                            

# Create an array that contains values between t0 and tEnd that have been 
# evenly spaced by dt
t=np.arange(t0,tEnd,dt)            

# Create an integer, n, that is the number of values in the array, t
n=int((tEnd-t0)/dt)                 

# Create arrays full of zeros that are n values long for each of the major variables
V=np.zeros(n)                      
N=np.zeros(n)                       
M=np.zeros(n)                       
H=np.zeros(n)                       
GNa=np.zeros(n)                     
GK=np.zeros(n)                      
INa=np.zeros(n)                     
IK=np.zeros(n)                      
IO=np.zeros(n)                      

# Set the initial value of the V, M, N, and H arrays to be equal to the initial conditions
V[0]=Vr                            
M[0]=M0                            
N[0]=N0                            
H[0]=H0                            

# Iterate through the t array --> Moving through time
for i in range(0,len(t)-1):
    
    # Calculate the variable GNa and the GK values at the specified time point
    GNa[i]=gNa*H[i]*M[i]**3     #Gna = (gNa max) * Hm^3
    GK[i]=gK*N[i]**4            #Gk = (gK max) * n^4   

    # Calculate the current
    # I = G * (Vm - equillibrium voltage)
    INa[i]=GNa[i]*(V[i]-ENa)    
    IK[i]=GK[i]*(V[i]-EK)       
    IO[i]=GO*(V[i]-EO)             
    
    # Set the relative vm value... 
    # This is used to find the rates of opening and closing the n gates
    vm=V[i]-Vr                     
     
    # Calculate the rate of opening/closing gates
    am=(0.1*(25-vm))/(np.exp((25-vm)/10)-1)     # rate of opening an activation gate
    bm=4*np.exp(-vm/18)                         # rate of closing an activation gate
    an=(0.01*(10-vm))/(np.exp((10-vm)/10)-1)    # rate of opening an n gate
    bn=0.125*np.exp(-vm/80)                     # rate of closing an n gate
    ah=0.07*np.exp(-vm/20)                      # rate of opening an inactivation gate
    bh=1/(np.exp((30-vm)/10)+1)                 # rate of closing an inactivation gate
    
    # Increment the next n, m, and h values by adding the old value to...
    # dn = dt*(am(1-m)-bm*m) (same formula for M and H)
    N[i+1]=N[i]+dt*(an*(1-N[i])-bn*N[i]); 
    M[i+1]=M[i]+dt*(am*(1-M[i])-bm*M[i]);            
    H[i+1]=H[i]+dt*(ah*(1-H[i])-bh*H[i]);            
    
    # if t is between 7 and 7.2
    # Solve for the next electric potential but add a stimulus current, Istim
    # Otherwise, proceed as normal
    if t[i]>=7 and t[i]<=7.2:                              
        V[i+1]=V[i]+dt*(-INa[i]-IK[i]-IO[i]+Istim)/Cm    
    else:
        V[i+1]=V[i]+dt*(-INa[i]-IK[i]-IO[i])/Cm         
  
    
  