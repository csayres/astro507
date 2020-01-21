import matplotlib.pyplot as plt
import numpy

n = 100 / 100*2 # particles per square cm
m = 1 # gram
k_b = 1.380648e-16 # erg/K
v = 100 / 5
vmin = 0
vmax = v * 3
ke = 0.5 * m * v**2
T = m * v / k_b

vs = numpy.linspace(vmin, vmax, 1000)
for tfact in [1,2,4,6,10,30,100]:
    TT = T*tfact
    p = n * m * vs / (2*numpy.pi*k_b*TT)*numpy.exp(-1*m*vs**2/(2*k_b*TT))
    plt.plot(vs, p, label="T_fact = %i"%tfact)
plt.legend()
plt.show()