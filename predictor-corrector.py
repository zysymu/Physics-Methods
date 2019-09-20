import matplotlib.pyplot as plt

global tf, x_in
tf = 50000
deltas = [0.1, 1, 10, 100, 1000, 10000]
x_in = 0.01
#general differential equation that is getting integrated:
#func := dx/dt = 0.0001*x*(2-x)

plt.style.use("ggplot")

#defining functions:
def f(x):
    return 0.0001*x*(2-x)


#plotting function
def plot(step, plot_x, plot_y, met): #plot_x e plot_y are lists
    for i in range(len(step)):
        plt.plot(plot_x[i], plot_y[i], label=str(step[i]) +" "+ met)


#method:
#if predCorr == False, only the first order (s=1) Adams-Bashforth method is used
#if predCorr == True, then the second order (s=2) Adams-Moulton method is implemented in order to make corrections
def adamsBash(dt, x_in, predCorr):
    ind = []
    x1 = x_in
    x0 = x_in
    t = 0
    while (t<=tf):
        x2 = x1 + dt*((3/2)*f(x1) - (1/2)*f(x0))
        
        if (predCorr==True):
            x2 = x1 + dt*((5/12)*f(x2) + (2/3)*f(x1) - (1/12)*f(x0))
        
        x0 = x1
        x1 = x2

        ind.append(x2)

        t += dt

    return ind



def time(dt, t):
    tind = []
    while (t<=tf):
        tind.append(t)
        t = t+dt
    
    return tind


#making a list of the values for different step sizes
ab_graph = [adamsBash(dt, x_in, predCorr=False) for dt in deltas]
pred_corr = [adamsBash(dt, x_in, predCorr=True) for dt in deltas]
t = [time(dt, t=0) for dt in deltas]


plot(deltas, t, ab_graph, met="")
plt.legend()
plt.title("Only predictor")
plt.show()


plot(deltas, t, pred_corr, met="")
plt.legend()
plt.title("Predictor-Corrector")
plt.show()

