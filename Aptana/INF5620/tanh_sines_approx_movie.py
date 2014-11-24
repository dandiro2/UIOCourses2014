import numpy as np
from math import pi, tanh, sin
import matplotlib.pyplot as plt

from scipy.integrate import quad
import sympy as sp
import matplotlib.animation as animation
fig, ax = plt.subplots()


def approximation(N):
    """
    Approximate exp(-x) by 
    1) using Galerkin or projection method f. The projection method 
    gives the same result as the least square method.
    """

    A = np.zeros((N+1,N+1)) 
    b = np.zeros(N+1)
     

    t = sp.symbols("t")      
 
    for i in range(N+1):
        phiiAndf = lambda x: sin( (2*i+1)*x )*tanh( 20*(x-pi)  )
        b[i],err =quad( phiiAndf,0,2*pi  ) 
         
        for j in range(N+1):
            phiiAndphij = lambda x: sin( (2*j+1)*x )*sin( (2*i+1)*x )
            #b[i] = sp.integrate(phi[i]*f, (x,0,2*pi) )
            
            A[i,j],err = quad(phiiAndphij, 0,2*pi  )
 
    #solve Ac = b for c
    c = np.linalg.solve(A, b)
       
    # approximation of exp(-t) by projection method
    phi = [sp.sin( (2*i+1)*t ) for i in range(N+1)]
    u = sum([ c[i]*phi[i] for i in range(N+1)])
    h = sp.lambdify(t,u,"numpy")
    z = np.linspace(0,2*pi,100)
    return h(z)
    
def functionTanh(s):
    z = np.linspace(0,2*pi,100)
    y = np.tanh(s*(z-pi))
    return y

# global variable pause to stop animation when N is to large
# when N is too large the integration diverge 
pause = False
 
def anime():
    """
    Animation function when N get larger
    """
    x = np.linspace(0, 2*np.pi, 100)        # x-array
    line, = ax.plot(x, np.tanh((x-pi)) )
 
 
    def animate(i):
        """
        matplotlib animation function
        """
        global pause
        if not pause:
            #line.set_ydata(functionTanh(i)   )  
            line.set_ydata(approximation(i))
        # stop if N > 23
        if i > 23:
            pause = True
        return line,
 
    def init():
        line.set_ydata(np.ma.array(x, mask=True))
        return line,
 
    ani = animation.FuncAnimation(fig, animate, init_func=init, frames = 1000,
        interval=10, blit=False, repeat = True)
    plt.xlim([0,2*pi])
    plt.ylim([-2,2])

    #y = np.tanh(2*(x-pi))
 
    #plot exp(-x) along with it approximation as N get larger
    #plt.plot(x,y,"r")
 
    plt.show()

if __name__=="__main__":
      
    anime()
    
    
    
    
    
    
    
    