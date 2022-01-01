# 7 FD methods for solving the 1d advection equation on the interval [0,LL]
# run the code. Type "python advection1d.py"
# The results can be visualized using gnuplt.
# On Linux system, type
#     gnuplot
#     load 'gpm'

import numpy as np
import math
import time
import sys



nMax = 10000
nbc = 2


print ("\n  7 Finite Difference methods for solving the 1d advection equation on the interval [0,LL]")
print ("  advection equation du/dt + du/dx = 0 \n") 
print ("  enter a number between 1 and 7 to choose a method \n") 
choice = int(input("upwind = 1, central scheme = 2, lax-wendrof = 3, beam-warming = 4, Fromm = 5, superbee = 6, minimod = 7: "))

if choice == 1:
    # (1) upwind
    numeric = "Upwind"
elif choice == 2:
    # (2)central scheme
    numeric = "Central"
elif choice == 3:
    # (3) lax-wendrof
    numeric = "Lax-Wendroff"
elif choice == 4:
    # (4) beam-warming
    numeric = "Beam-Warming"
elif choice == 5:
    # (5) fromms
    numeric = "Fromm"
elif choice == 6:
    # (6) superbee
    numeric = "Superbee"
    numeric = "TVD"
elif choice == 7:
    # (7) minimod
    numeric = "Minimod"

#nn = input("read nn: ")
nn = 50
# interval [0, LL]
LL = 10.0
cfl = 0.5
dx = LL/nn
dt = dx*cfl
adtdx = dt/dx
dtdx = dt/dx

print("dx = %g, dt %g\n" %(dx, dt))
print("adtdx = %g, cfl = %g\n" %(adtdx, cfl))

tmax = int(4.0/dt)
print("tmax = %g, time %g\n" %(tmax, tmax*dt))

#initial condition
#phi = []
#fl = []
phi = np.zeros(10000)
phi = np.zeros(nn+2*nbc+1)
fl = np.zeros(10000)
L0 = 1.0
#for i in range(-nbc, nn+nbc+1):
#print("%d %d\n" %(nbc, nn))
for i in range(0, nn+2*nbc+1):
    xx = (i-nbc+0.5)*dx
    #print("%d" %i)
    if xx <= 1.0:
        phi[i] = 1.0
    else:
        phi[i] = 0.0

#time marching forward
gpMovie =  open("gpm", "w")
ims = []

def minimod(a,b):
    r = a*b
    if r <=0:
        return 0
    elif abs(a) < abs(b):
        return a
    else:
        return b

def maxmod(a,b):
    r = a*b
    if r<=0:
        return 0
    if abs(a) > abs(b):
        return a
    else:
        return b

# fl is the flux in the FD scheme. !!!

# to complete this program, define the flux 
# for (5) FROMM and (6) Superbee
# You can immitate (3) ax-wendrof and (7) minimod
# or Write your own formulas

for it in range(1, tmax+1):
    for i in range(nbc, nn+nbc+1):
        if choice == 1:
            # (1) upwind
            fl[i] = phi[i-1]
        elif choice == 2:
            # (2)central scheme 
            fl[i] = 0.5*(phi[i-1]+phi[i])
        elif choice == 3:
            # (3) lax-wendrof
            #fl[i] = 0.5*(phi[i-1]+phi[i]) + 0.5*adtdx*(phi[i-1]-phi[i])
            fl[i] = phi[i-1] + 0.5*(1.0-adtdx)*(phi[i]-phi[i-1])
        elif choice == 4:
            # (4) beam-warming
            fl[i] = phi[i-1] + 0.5*(1.0-adtdx)*(phi[i-1]-phi[i-2])
        elif choice == 5:
            # (5) fromms
            fl[i] = phi[i-1] + 0.25*(1.0-adtdx)*(phi[i]-phi[i-2])
        elif choice == 6:
            # (6) superbee
            #print ("Complete the Superbee scheme. exit\n")
            #sys.exit() 
            lmter = maxmod(minimod(phi[i-1]-phi[i-2], 2*(phi[i]-phi[i-1]))  ,   minimod(2*(phi[i-1]-phi[i-2]), (phi[i]-phi[i-1])))
           # lmter = max( min(phi[i-1]-phi[i-2], 2*(phi[i]-phi[i-1]))  ,   min(2*(phi[i-1]-phi[i-2]), (phi[i]-phi[i-1])))
            fl[i] = phi[i-1] + 0.5*(1.0-adtdx)*lmter
        elif choice == 7:
            # (7) minimod
            lmter = minimod(phi[i]-phi[i-1], phi[i-1]-phi[i-2])
            fl[i] = phi[i-1] + 0.5*(1.0-adtdx)*lmter

#  update phi
    for i in range(nbc, nn+nbc):
       # phi[i] = phi[i]+adtdx*(fl[i]-fl[i+1])
       phi[i] = phi[i]+dtdx*(fl[i]-fl[i+1])

    if it < 10:
        filename = numeric +".0"+ str(it)
    elif it < 100:
        filename = numeric +"."+ str(it)
        
    if it < 10:
        exact_name = "exact.0"+str(it) 
    elif it< 100:
        exact_name = "exact."+str(it) 

    #print(filename+"\n")
    #exit()

    file = open(filename, "w")
    x = np.zeros(10000)
    xval = np.zeros(nn+2*nbc+1)
    for i in range(0, nn+2*nbc+1):
        xx = (i-nbc+0.5)*dx
        xval[i] = (i-nbc+0.5)*dx
        file.write("%g %g\n" %(xx, phi[i]))
    #file.close()

    gpMovie.write("plot [0:10][-.5:1.5]'"+filename+"'w p pt 5 ps 2,'"+filename+"'w l lw 2,'"+exact_name+"'w l lw 3\n")
    #gpMovie.write("plot [0:10][-.5:1.5]'"+numeric+"'w p pt 7 ps 2,'exact.%d'w l lw 3\n" %(it))
    #gpMovie.write("plot [0:10][-.5:1.5]'superbee.%d'w p pt 7 ps 2,'exact.%d'w l lw 3\n" %(it, it))

    xex = np.zeros(4)
    yex = np.zeros(4)
    xex[0] = 0
    yex[0] = 1
    xex[1] = L0 + dt*it
    yex[1] = 1
    xex[2] = L0 + dt*it
    yex[2] = 0
    xex[3] = LL
    yex[3] = 0


    file = open(exact_name, "w")
    file.write("0 1\n")
    file.write("%g 1\n" %(L0 + dt*it))
    file.write("%g 0\n" %(L0 + dt*it))
    file.write("%g 0\n" %LL)
    file.close()


