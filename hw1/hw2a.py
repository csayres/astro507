import matplotlib.pyplot as plt
import numpy
from cParticle import Box
import pickle

k_b = 1.380648e-16 # erg/K
width = 100*2
height = 100*2
radius = 2
mass = 1
nParticles = 400*4
v = width / 10
ke = 0.5 * mass * v**2
T = ke / k_b


def runSim():
    seed = 0
    box = Box(width, height, seed)
    for ii in range(nParticles):
        box.addRandomParticle(mass, radius, v)
    steps = 10000
    dt = 0.02 * radius / v
    saveEvery = 1
    box.runSim(steps, dt, saveEvery)
    # data = numpy.array(box.particleSteps)
    pickle.dump(box.particleSteps, open("ds.p", "wb"))


# runSim()
data = pickle.load(open("ds.p", "rb"))
bins = numpy.linspace(0,3*v,50)
histArr = []
relaxed = numpy.asarray(data[600:])
speeds = relaxed[:,:,4].squeeze()
# for speed in relaxed:
#     histArr.append(numpy.histogram(speed, bins)[0])
# histArr = numpy.asarray(histArr)
# histAvg = numpy.mean(histArr, axis=0)
plt.figure(figsize=(10,10))
plt.hist(speeds.flatten(), bins=bins, density=True, alpha=0.7)
vs = numpy.linspace(0, 3*v, 10000)
p = mass * vs / (k_b * T) * numpy.exp(-1 * mass * vs**2 / (2 * k_b * T))
dv = numpy.diff(vs)[0]
pSum = numpy.sum(p) * dv
plt.plot(vs, p/pSum, '--r')
plt.xlabel("velocity (cm/s)")
plt.savefig("hw1a.png", dpi=250)

