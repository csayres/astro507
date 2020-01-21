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
            [particle1.xPos - particle2.xPos, particle1.yPos - particle2.yPos]
        )
        return dist < (particle1.radius + particle2.radius)

    def handleCollision(self, particle1, particle2):
        """Modify the velocities of particle 1 and 2 for an elastic collision

        https://en.wikipedia.org/wiki/Elastic_collision
        """
        contactAngle = numpy.arctan2(
            particle1.yPos - particle2.yPos,
            particle1.xPos - particle2.xPos
        )
        A1 = (
            particle1.speed *
            numpy.cos(particle1.motionAngle - contactAngle) *
            (particle1.mass - particle2.mass) +
            2 * particle2.mass * particle2.speed *
            numpy.cos(particle2.motionAngle - contactAngle)
        ) / (particle1.mass + particle2.mass)

        A2 = (
            particle2.speed *
            numpy.cos(particle2.motionAngle - contactAngle) *
            (particle2.mass - particle1.mass) +
            2 * particle1.mass * particle1.speed *
            numpy.cos(particle1.motionAngle - contactAngle)
        ) / (particle2.mass + particle1.mass)

        # new x velocity for particle 1
        v1x = A1 * numpy.cos(contactAngle) + \
            particle1.speed * numpy.sin(particle1.motionAngle - contactAngle) * \
            numpy.cos(contactAngle + numpy.pi / 2)
        # new y velocity for particle 1
        v1y = A1 * numpy.sin(contactAngle) + \
            particle1.speed * numpy.sin(particle1.motionAngle - contactAngle) * \
            numpy.sin(contactAngle + numpy.pi / 2)

        # new x velocity for particle 2
        v2x = A2 * numpy.cos(contactAngle) + \
            particle2.speed * numpy.sin(particle2.motionAngle - contactAngle) * \
            numpy.cos(contactAngle + numpy.pi / 2)
        # new y velocity for particle 2
        v2y = A2 * numpy.sin(contactAngle) + \
            particle2.speed * numpy.sin(particle2.motionAngle - contactAngle) *\
            numpy.sin(contactAngle + numpy.pi / 2)

        particle1.xVel = v1x
        particle1.yVel = v1y
        particle2.xVel = v2x
        particle2.yVel = v2y

    def addParticle(self, particle):
        self.particles.append(particle)

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
            p = Particle(mass, radius, x, y, xVel, yVel, self.width, self.height)
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
                break # from while loop

        print("%i particles added"%len(self.particles))

    def runSim(self, maxSteps=10000, saveEvery=1):
        """Start the simulation, pick a reasonable timestep based on ke
        and smallest particle radius
        """
        self.timeSteps = []
        maxSpeed = 0
        minRadius = numpy.inf
        for particle in self.particles:
            if particle.speed > maxSpeed:
                maxSpeed = particle.speed
            if particle.radius < minRadius:
                minRadius = particle.radius

        # set the time step such that dt*maxSpeed < 0.01 the particle radius
        self.dt = 0.005 * minRadius / maxSpeed

        for step in range(maxSteps):
            print("step", step)
            # update positions
            for particle in self.particles:
                # updatePostion handles reversal at walls
                particle.updatePosition(self.dt)

            # look at all pairs of particles for collisions
            for particle1, particle2 in itertools.combinations(self.particles, 2):
                if self.isCollided(particle1, particle2):
                    self.handleCollision(particle1, particle2)

            # save state
            if step % saveEvery == 0:
                particleStates = []
                for particle in self.particles:
                    particleStates.append(
                        [particle.xPos, particle.yPos,
                         particle.mass, particle.speed]
                    )
                self.timeSteps.append(particleStates)








