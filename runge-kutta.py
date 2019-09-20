"""
Code authored by zysymu~
Made with love for Metodos Computacionais da Fisica B.
It contains different Runge-Kutta methods, and also a practical way to plot the results using different
step sizes in order to make graphical comparisons.
"""
import matplotlib.pyplot as plt
from math import *
import numpy as np

plt.style.use("ggplot")

#defining variables
global x_in, tau, t, tf
steps = [0.1, 0.5] #can contain different step sizes in order to make an automatic comparison
tf = 6
x_in = 10.
tau = 2.

#analytical function (it`s going to be used in order to make comparisons with the algorithm)
#also, it's a simple radioactive decay equation
ta = np.linspace(0, tf)
xa = x_in * np.exp(-ta/tau)


#general plotting function
def plot(step, plot_x, plot_y, met): #plot_x and plot_y are lists containing sublists
    for i in range(len(step)):
        plt.plot(plot_x[i], plot_y[i], label=str(step[i]) +" "+ met)
 

def time(dt, t):
    tind = []
    while (t<=tf):
        tind.append(t)
        t = t+dt
    
    return tind

#runge-kutta 2 (rk2) and runge-kutta 4 (rk4) algorithms:
def rk2_heun(dt, x):
    ind = []
    t=0
    while (t<=tf):
        k1 = -(x/tau)
        x_aux = x + k1*dt
        k2 = -(x_aux/tau)
        ind.append(x)
        x = x + (1/2)*(k1+k2)*dt
        t = t+dt
        
    return ind


def rk2_midpoint(dt, x):
    ind = []
    t=0
    while (t<=tf):
        k1 = -(x/tau)
        x_aux = x + k1*(dt/2)
        k2 = -(x_aux/tau)
        ind.append(x)
        x = x + k2*dt
        t = t+dt
        
    return ind


def rk2_ralston(dt, x):
    ind = []
    t=0
    while (t<=tf):
        k1 = -(x/tau)
        x_aux = x + (3/4)*k1*dt
        k2 = -(x_aux/tau)
        ind.append(x)
        x = x + (1/3)*k1*dt + (2/3)*k2*dt
        t = t+dt
        
    return ind



def rk4_3o8(dt, x):
    ind = []
    t=0
    while (t<=tf):
        ind.append(x)
        k1 = -(x/tau)
        k2 = -(x + k1*(dt/3.))/tau
        k3 = -(x - k1*(dt/3.) + k2*dt)/tau
        k4 = -(x + k1*dt - k2*dt + k3*dt)/tau
        x = x + (1./8.)*(k1 + 3.*k2 + 3.*k3 + k4)*dt
        t = t+dt
        
    return ind
        


def rk4_classic(dt, x):
    ind = []
    t=0
    while (t<=tf):
        ind.append(x)
        k1 = -(x/tau)
        k2 = -(x + (1./2.)*k1*dt)/tau
        k3 = -(x + (1./2.)*k2*dt)/tau
        k4 = -(x + k3*dt)/tau
        x = x + (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4)*dt
        t = t+dt
        
    return ind



#applying
#rk4 classic
rk4_plot = [rk4_classic(dt, x_in) for dt in steps]
    
#rk4 3/8
rk4_3o8 = [rk4_3o8(dt, x_in) for dt in steps]

#plots
t_plot = [time(dt, t=0) for dt in steps]    
midpoint_plot = [midpoint(dt, x_in) for dt in steps]

plt.plot(ta, xa, label="analytic", c="black")
plot(steps, t_plot, midpoint_plot, met="midpoint")
plot(steps, t_plot, rk4_plot, met="rk4")

plt.xlabel("time")
plt.ylabel("num of atoms")
plt.legend()
plt.show()

