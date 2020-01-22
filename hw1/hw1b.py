from particles import Box
import numpy
numpy.random.seed(7)
from shapely.geometry import Point
from descartes import PolygonPatch
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
from subprocess import Popen
import glob
import os


width = 100
height = 100
radius = 2
mass = 1
nParticles = 100
v = width / 5 # 5 seconds to cross box
ke = 0.5 * mass * v**2
box = Box(width, height)
for ii in range(nParticles):
    box.addRandomParticle(mass, radius, ke)
steps = 50000
box.runSim(steps, saveEvery=50)

def plotOne(step):
    fig = plt.figure(figsize=(6,8))
    ax1 = plt.subplot2grid((4,3), (0,0), rowspan=3, colspan=3)
    ax2 = plt.subplot2grid((4,3), (3,0), colspan=3)
    for ii, particle in enumerate(box.simSteps[step]):

        if ii in box.collisionFlags[step]:
            topcolor="red"
        else:
            topcolor="blue"
        pt = Point(particle[2], particle[3]).buffer(particle[1], cap_style=1)
        patch = PolygonPatch(pt, fc=topcolor)
        ax1.add_patch(patch)
    ax1.set_ylim([0, height])
    ax1.set_xlim([0, width])

    # plot histogram of speeds
    speeds = numpy.linalg.norm(box.simSteps[step][:,4:6], axis=1)
    bins = numpy.linspace(0,3*v,25)
    ax2.hist(speeds, bins=bins)
    ax2.set_xlim([0, 3*v])
    # ax2.set_ylim([0, nParticles*1.1])
    ax2.set_xlabel("velocity")
    plt.savefig("step_%08d.png"%step, dpi=150)
    plt.close()

# plotOne(0)


p = Pool(cpu_count())
p.map(plotOne, range(len(box.simSteps)))

fps = 60
args = ['ffmpeg', '-r', '%i'%fps, '-f', 'image2', '-i', 'step_%08d.png',
    '-pix_fmt', 'yuv420p', 'moviehw1b.mp4']
movie = Popen(args)
movie.wait()
# clean up imgs
imgs = glob.glob("step*.png")
for img in imgs:
    os.remove(img)
