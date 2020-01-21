from particles import Box
import numpy
numpy.random.seed(0)
from shapely.geometry import Point
from descartes import PolygonPatch
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
from subprocess import Popen
import glob
import os


width = 100
height = 100
radius = 10
mass = 1
v = width / 5 # 5 seconds to cross box
ke = 0.5 * mass * v**2
box = Box(width, height)
box.addRandomParticle(mass, radius, ke)
box.addRandomParticle(mass*10, radius, ke)
steps = 100000
box.runSim(steps, saveEvery=100)

def plotOne(step):
    plt.figure(figsize=(10,10))
    ax = plt.gca()
    for particle in box.timeSteps[step]:
        if particle[2] < 2:
            topcolor="blue"
        else:
            topcolor="red"
        pt = Point(particle[0], particle[1]).buffer(radius, cap_style=1)
        patch = PolygonPatch(pt, fc=topcolor)
        ax.add_patch(patch)
    ax.set_ylim([0, height])
    ax.set_xlim([0, width])
    plt.savefig("step_%08d.png"%step, dpi=150)
    plt.close()

# plotOne(1)

p = Pool(cpu_count())
p.map(plotOne, range(len(box.timeSteps)))

fps = 50
args = ['ffmpeg', '-r', '%i'%fps, '-f', 'image2', '-i', 'step_%08d.png',
    '-pix_fmt', 'yuv420p', 'twoParticles.mp4']
movie = Popen(args)
movie.wait()
# clean up imgs
imgs = glob.glob("step*.png")
for img in imgs:
    os.remove(img)