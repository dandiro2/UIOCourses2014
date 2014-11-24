


import sympy as sp
import numpy as np
from plot_test import *

def test_constant():
    I = lambda x,y: 0.3  # constant initial elevation
    V = lambda x,y: 0  # du/dt at t=0
    q = lambda x,y: 1
    f = lambda x,y,t:0
    #dx = 0.05; dy = 0.05; dt = 0.05
    #solve(I,q,f,V,b,dt,dx,dy,Lx,Ly,T)
    dt=0.01;dx=0.05;dy=0.05;T=1.0;Lx=1.0;Ly=1.0
    N = 21
    u_num = solve(I,q,f,V,0,dt,dx,dy,Lx,Ly,T)  
    sol = 0.3*np.ones((N,N))    # should be 0.3 in all points
    tol = 1e-13
    success = abs(np.max(u_num-sol) < tol)

    
    print 'max error constant: %e ' %np.max(u_num-sol)
    assert success



def test_plug():
    # the plug should return to its initial state after this simulation
    def Ix(x,y):
        a = (x<0.75)
        b = (x>0.25)    #booleans
        return 0.3*a*b
    def Iy(x,y):
        a = (y>0.4)
        b = (y<0.6)
        return 0.3*a*b    # plug in y direction
     
    V = lambda x,y: 0
    q = lambda x,y: 1
    f = lambda x,y,t: 0
    b = 0
    dx = 0.05; dy = 0.05; dt = 0.05
    plug_x = solve(Ix,f,q,V,b)
    plug_y = solve(Iy,f,q,V,b)
     
    error_x = np.max(abs(plug_x.u[-1]-plug_x.u[0]))
    error_y = np.max(abs(plug_y.u[-1]-plug_y.u[0]))
    tol = 1e-13
    success = (error_x < tol and error_y < tol)
    if __name__=="__main__":
        plug_y.animate(filename = 'plug_y.gif')
    print 'max error plug: %e' %(max(error_x,error_y))
    assert success
 
def test_standing():
    T = 1.0; L_x = 2.; L_y = 1.    
    k_x = 4*np.pi/L_x; k_y = 2*np.pi/L_y
 
    A = 0.3
    I = lambda x,y : A*np.cos(k_y*y)*np.cos(k_x*x)
    q = lambda x,y : 1
    V = lambda x,y : 0
    f = lambda x,y,t: 0
     
    b = 0
    wave = wavesolver(I,f,q,V,b)    # create wavesolver object
    N_test = 5
    H = [0.2*(0.5**i) for i in range(N_test)]
    E = []
    for i in range(N_test):
        h = H[i]
        dt = 0.5*h; dx = h; dy = h
         
        M = int(round(T/dt)) + 1; N_x = int(round(L_x/dx)) + 1; N_y = int(round(L_y/dy)) + 1
 
        wave.solve(dt,dx,dy,T,L_x,L_y)    # solve with given h
        omega = np.sqrt(k_x**2 + k_y**2)
        x = np.linspace(0,L_x,N_x)
        y = np.linspace(0,L_y,N_y)
        t = np.linspace(0,T,M)
        X,Y = np.meshgrid(x,y)
        X = X.transpose()
        Y = Y.transpose()    # gives correct shapes
        exact = np.zeros((M,N_x,N_y))
        for n in range(M):
            exact[n,:,:] = (A*np.cos(k_x*X)*np.cos(k_y*Y))*np.cos(omega*t[n])    # exact solution
 
        error = abs(exact - wave.u)
        max_error = np.max(error)
        E.append(max_error)
        if i>0:
            r = np.log(E[i]/E[i-1])/np.log(H[i]/H[i-1])
        if __name__=="__main__" and i == 1:   # example of standing solution
                wave.animate(filename = 'standing.gif')
    success = (abs(r-2) < 0.01)
    print 'standing wave convergence: r_%d = %f ' %(i,r)
    assert success
 
def test_manufracted():
    T = 2.; L_x = 2.; L_y = 1.
    k_x = 5*np.pi/L_x; k_y = 1*np.pi/L_y
     
 
    # Choose arbitrary parameters/functions
    A = 0.3; B = 0.3; b = 0.01;c = b/2.
    x_s,y_s,t_s = sp.symbols('x_s,y_s,t_s')
    q_s = 0.3*x_s*y_s
    omega = sp.sqrt(k_x**2+k_y**2)
 
    # exact solution    
    u_e = (A*sp.cos(omega*t_s) + B*sp.sin(omega*t_s))*sp.exp(-c*t_s)*(sp.cos(k_x*x_s)*sp.cos(k_y*y_s))
     
    #sympy calculations for f,I and V
    u_tt = sp.diff(sp.diff(u_e,t_s),t_s)
    u_t = sp.diff(u_e,t_s)
    qu_x_x = sp.diff(q_s*sp.diff(u_e,x_s),x_s)
    qu_y_y = sp.diff(q_s*sp.diff(u_e,y_s),y_s)
 
    f_s = u_tt + b*u_t - qu_x_x - qu_y_y
    V_s = sp.diff(u_e,t_s).subs(t_s,0)
    I_s = u_e.subs(t_s,0)
 
    #turning expressions to functions
    u_func = sp.lambdify([x_s,y_s,t_s],u_e,modules='numpy')
    q = sp.lambdify([x_s,y_s],q_s,modules='numpy')
    I = sp.lambdify([x_s,y_s],I_s,modules='numpy')
    V = sp.lambdify([x_s,y_s],V_s,modules='numpy')
    f = sp.lambdify([x_s,y_s,t_s],f_s,modules='numpy')
 
    # make a wavesolver object
    wave = wavesolver(I,f,q,V,b)
    N_test = 5
    H = [0.2*(0.5**i) for i in range(N_test)]
    E = []
    for i in range(N_test):
        h = H[i]
        dt = 0.5*h; dx = h; dy = h        # 0.5 on dt for safety    
         
        M = int(round(T/dt)) + 1; N_x = int(round(L_x/dx)) + 1; N_y = int(round(L_y/dy)) + 1
        wave.solve(dt,dx,dy,T,L_x,L_y) # solve with given h
 
        t = np.linspace(0,T,M)
        x = np.linspace(0,L_x,N_x);y = np.linspace(0,L_y,N_y)
        X,Y = np.meshgrid(x,y)
        X = X.transpose()
        Y = Y.transpose()
        exact = np.zeros((M,N_x,N_y))
     
        for n in range(M):
            exact[n,:,:] = u_func(X,Y,t[n])
         
     
        error = np.max(abs(exact-wave.u))
        E.append(error)
        if i>0:
            r = np.log(E[i]/E[i-1])/np.log(H[i]/H[i-1])
        if __name__=="__main__" and i == 1:
                wave.animate(filename='manufracted.gif')
    # satisfied with abs(r-2)<0.01 
    success = (abs(r-2) < 0.01)
    print 'manufracted solution convergence: r_%d = %f ' %(i,r)
    assert success
             

if __name__=='__main__':
    test_constant()
    #test_plug()
    #test_standing()
    #test_manufracted()
