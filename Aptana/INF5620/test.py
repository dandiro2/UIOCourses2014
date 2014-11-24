"""
Animating the approximation
"""
from tanh_sines_approx1 import approximation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import pi
fig, ax = plt.subplots()

x = np.linspace(0, 2*np.pi, 100)        # x-array
line, = ax.plot(x, np.tanh(20*(x-pi)) )

pause = False
def animate(i):
    global pause
    if not pause:
        line.set_ydata(approximation(i)   )  
    # stop if N > 23
    if i > 23:
        pause = True
    return line,

def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,

ani = animation.FuncAnimation(fig, animate, init_func=init, frames = 1000,
    interval=900, blit=False, repeat = True)
plt.xlim([0,2*pi])
plt.ylim([-2,2])
y = np.tanh(20*(x-pi))

plt.plot(x,y,"r")

plt.show()
