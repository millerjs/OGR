from pylab import *
from numpy import *
from scipy import *
from scipy import optimize
import math


x,y  = loadtxt('out.txt',unpack=True, usecols=[0,1], skiprows=5)
plt = matplotlib.pyplot.figure()
ax = axes()
ax.plot(x,y,"bo")
ax.set_title("Reproduced Plot")
plt.set_facecolor('white')
show()
