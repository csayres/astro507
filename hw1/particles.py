import itertools
import numpy

class Particle(object):
    def __init__(self, mass, radius, xPos, yPos, xVel, yVel, maxX, maxY):
        """
           mass in g
           radius in cm
           xPos, yPos cm
           xVel, yVel cm/s
        """
        self.mass = mass
        self.radius = radius
        self.xPos = xPos
        self.yPos = yPos
        self.xVel = xVel
        self.yVel = yVel
        self.maxX = maxX
        self.maxY = maxY

    def updatePosition(self, dt):
        """Propagate particle according to velocities
        reverse velocity if we're at max (wall boundary)
        """
        self.xPos = self.xPos + self.xVel * dt
        self.yPos = self.yPos + self.yVel * dt

        # reverse direction if we're at a wall
        if self.xPos + self.radius > self.maxX and self.xVel > 0:
            self.xVel = -1 * self.xVel
        if self.xPos - self.radius < 0 and self.xVel < 0:
            self.xVel = -1 * self.xVel
        if self.yPos + self.radius > self.maxY and self.yVel > 0:
            self.yVel = -1 * self.yVel
        if self.yPos - self.radius < 0 and self.yVel < 0:
            self.yVel = -1 * self.yVel

    @property
    def speed(self):
        """Speed of particle"""
        return numpy.sqrt(self.xVel**2 + self.yVel**2)

    @property
    def p(self):
        """momentum of particle in x and y"""
        return self.mass * self.xVel, self.mass * self.yVel

    @property
    def ke(self):
        """kinetic energy of particle"""
        return 0.5 * self.mass * (self.xVel**2 + self.yVel**2)

    @property
    def motionAngle(self):
        """Return the direction of motion as an angle
        """
        return numpy.arctan2(self.yVel, self.xVel)


class Box(object):
    def __init__(self, width, height):
        """input width / height in cm
        """
        # 0,0 is box center
        self.width = width
        self.height = height
        self.particles = []
        self.dt = None # set at simulation time

    def randomPoint(self, radius):
        """Return a random xy point inside of box
        """
        x = numpy.random.uniform(0+radius, self.width-radius)
        y = numpy.random.uniform(0+radius, self.height-radius)
        return x, y

    def isCollided(self, particle1, particle2):
        """Return True if these particles are collided
        """
        dist = numpy.linalg.norm(
            [particle1[2] - particle2[2], particle1[3] - particle2[3]]
        )
        return dist < (particle1[1] + particle2[1])


    def findCollision(self, pIndex):
        """Check if this particle is collided, expect that is
        may only be collided with 1 other particle

        return the set {pIndex, collidedNeighbor}
        """

        # calculate distance to all particles
        xys = self.particles[:,2:4]
        radii = self.particles[:,1]
        dxys = xys - xys[pIndex]
        dist = numpy.linalg.norm(dxys,axis=1)
        # note we can't handle differing radii
        collidedIndices = set(numpy.argwhere(dist <= 2*radii).flatten())
        nCollided = len(collidedIndices)
        if nCollided > 2:
            print("warning: > 2 body collision, shouldn't happen: %i"%nCollided)
        if nCollided == 1:
            # ignore a self-collision!
            return None
        # returns, {pIndex, collidedIndex}
        return set(collidedIndices)


    def handleCollision(self, p1index, p2index):
        """Modify the velocities of particle 1 and 2 for an elastic collision

        arguments are indices in particle list

        https://en.wikipedia.org/wiki/Elastic_collision
        """


        # contact angle
        ca = numpy.arctan2(
            self.particles[p1index][3] - self.particles[p2index][3],
            self.particles[p1index][2] - self.particles[p2index][2]
        )

        # calculate motion angles for particles
        ma1 = numpy.arctan2(self.particles[p1index][5], self.particles[p1index][4])
        ma2 = numpy.arctan2(self.particles[p2index][5], self.particles[p2index][4])

        # calculate speed for particles
        sp1 = numpy.linalg.norm(self.particles[p1index][4:6])
        sp2 = numpy.linalg.norm(self.particles[p2index][4:6])

        # masses of particles
        m1 = self.particles[p1index][0]
        m2 = self.particles[p2index][0]


        A1 = (
            sp1 * numpy.cos(ma1 - ca) * (m1 - m2) +
            2 * m2 * sp2 * numpy.cos(ma2 - ca)
        ) / (m1 + m2)

        A2 = (
            sp2 * numpy.cos(ma2 - ca) * (m2 - m1) +
            2 * m1 * sp1 * numpy.cos(ma1 - ca)
        ) / (m2 + m1)

        # new x velocity for particle 1
        v1x = A1 * numpy.cos(ca) + sp1 * numpy.sin(ma1 - ca) * numpy.cos(ca + numpy.pi / 2)
        # new y velocity for particle 1
        v1y = A1 * numpy.sin(ca) + sp1 * numpy.sin(ma1 - ca) * numpy.sin(ca + numpy.pi / 2)

        # new x velocity for particle 2
        v2x = A2 * numpy.cos(ca) + sp2 * numpy.sin(ma2 - ca) * numpy.cos(ca + numpy.pi / 2)
        # new y velocity for particle 2
        v2y = A2 * numpy.sin(ca) + sp2 * numpy.sin(ma2 - ca) * numpy.sin(ca + numpy.pi / 2)

        self.particles[p1index][4] = v1x
        self.particles[p1index][5] = v1y
        self.particles[p2index][4] = v2x
        self.particles[p2index][5] = v2y

    def addRandomParticle(self, mass, radius, ke):
        """Add a random particle to the box, give it a random position
        and random direction with kinetic energy specified
        id : an identifier for this particle
        mass: mass of particle in grams
        radius: radius of particle in cm
        ke: kinetic energy in gcm2/s2 (erg)
        """
        speed = numpy.sqrt(2 * ke / mass)
        # pick a random angle 0 - 2pi
        ang = numpy.random.uniform(0, 2*numpy.pi)
        xVel = numpy.cos(ang) * speed
        yVel = numpy.sin(ang) * speed
        # find a non-colliding place for it in the box
        while True:
            x, y = self.randomPoint(radius)
            p = [mass, radius, x, y, xVel, yVel]
            if not self.particles:
                print('added first particle')
                # this is first particle
                self.particles.append(p)
                break # from while loop

            # check if this particle collides with any others
            isCollided = False
            for otherParticle in self.particles:
                if self.isCollided(p, otherParticle):
                    isCollided = True
                    break

            if not isCollided:
                # particle isn't collided, save it
                self.particles.append(p)
                print('added %i particle'%(len(self.particles)))
                break # from while loop


    def updatePositions(self):
        """Propagate all particles according to their velocities,
        dt is specified in runSim()
        """
        self.particles[:, 2] = self.particles[:, 2] + self.particles[:, 4] * self.dt
        self.particles[:, 3] = self.particles[:, 3] + self.particles[:, 5] * self.dt

        # find particles at box edges and reverse their velocities
        # may want to check that velocities haven't been updated already?
        rightEdge = self.particles[:, 2] + self.particles[:, 1] >= self.width
        leftEdge = self.particles[:, 2] - self.particles[:, 1] <= 0
        topEdge = self.particles[:, 3] + self.particles[:, 1] >= self.height
        bottomEdge = self.particles[:, 3] - self.particles[:, 1] <= 0

        # modify x velocities for left and right edges
        self.particles[rightEdge, 4] = -1*self.particles[rightEdge, 4]
        self.particles[leftEdge, 4] = -1*self.particles[leftEdge, 4]

        # modify y velocities for left and right edges
        self.particles[topEdge, 5] = -1*self.particles[topEdge, 5]
        self.particles[bottomEdge, 5] = -1*self.particles[bottomEdge, 5]


    def runSim(self, maxSteps=10000, saveEvery=1, dt=None):
        """Start the simulation, pick a reasonable timestep based on ke
        and smallest particle radius
        """
        #
        self.simSteps = [] # holds snapshots of the simulation
        self.collisionFlags = []

        # convert particle list to numpy array for quicker computations
        self.particles = numpy.array(self.particles)

        if dt is None:
            # set the time step such that dt*maxSpeed < 0.01 the particle radius
            maxSpeed = numpy.max(numpy.linalg.norm(self.particles[:,4:6], axis=1))
            minRadius = numpy.min(self.particles[:,1])
            self.dt = 0.0005 * minRadius / maxSpeed
        else:
            self.dt = dt

        for step in range(maxSteps):
            print("step", step)
            # update positions propagates particles
            # it also reflects velocities at edges of box
            self.updatePositions()


            collidedPairs = []
            for ii in range(len(self.particles)):
                # findCollision returns a set
                collidedPair = self.findCollision(ii)
                if collidedPair is not None:
                    collidedPairs.append(collidedPair)

            # only operate on unique collisions
            # (we find both A-B and B-A collisions)
            collisionFlags = []
            for collidedPair in numpy.unique(collidedPairs):
                index1, index2 = list(collidedPair)
                self.handleCollision(index1, index2)
                collisionFlags.append(index1)
                collisionFlags.append(index2)


            # save state
            if step % saveEvery == 0:
                self.simSteps.append(numpy.copy(self.particles)) # ensures a copy
                self.collisionFlags.append(collisionFlags)




