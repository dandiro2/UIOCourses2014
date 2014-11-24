
import numpy as np



def solver(I,N,x,y):
    u = np.zeros((N,N))
    
    for i in range(N):
        for j in range(N):
            u[i,j] = I(x[i],x[j])
            
    print u
    

def test():
    N, L = 10,5
    U = lambda x,y: 0 if abs(x) < 5 else 10
    x = np.linspace(0,L,N)
    y = np.linspace(0,L,N)
    I = U
    
    V = solver(I,N,x,y)

    print V


test()
    




