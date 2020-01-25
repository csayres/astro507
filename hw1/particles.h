# pragma once
#include <vector>
#include <array>

class Particle{
public:
    double mass;  // g
    double radius; // cm
    double x; // cm
    double y; // cm
    double vx; // cm/s
    double vy; // cm/s
    Particle(double mass, double radius, double x, double y, double vx, double vy);
    // input mass, radius, x, y, velociy_x, velocity_y in cgs units
    void updatePosition(double dt);
    // propatage motion by dt (seconds)
    std::array<double, 2> nextPosition(double dt);
    // return [x,y] array for next postion given dt (seconds)
    void setVelocities(double vx, double vy);
    // update the particles velocities in x and y
    double speed();
    // return the particles speed in cm/s
    double direction();
    // return the angle of motion in radians
};

class Box{
public:
    double width;  // width of box cm
    double height; // height of box cm
    double dt; // time step for simulation
    std::vector<std::shared_ptr<Particle>> particles; // list of all particles in the box
    // list of saved states during the simulation
    std::vector<std::vector<std::array<double, 5>>> particleSteps;
    std::vector<double> pressureSteps; // at each step
    std::vector<double> collisionSteps; // collisions registered at each step
    Box(double width, double height, int seed=0);
    // input box width and height in cm, seed is the random seed
    std::array<double, 2> randomPoint(double radius);
    // return and [x,y] random point in the box accounting for particle radius (cm)
    void addRandomParticle(double mass, double radius, double speed);
    // add a random particle with given mass radius and speed
    // ensure it doesn't collide upon placement.
    bool isCollided(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);
    // return true if particles are collided, and approacing eachother false otherwise
    void handleCollision(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);
    // input pInd1, pInd2 are the indices of particles in particles.  Update their
    // velocities for the case of an elasitc collision
    void runSim(int steps, double dt, int saveEvery);
    // Begin the simulation, specify the number of steps and time step dt (seconds),
    // saveEvery specifies the frequency (in steps) at which to save the system
    // state into the particleSteps array.
    std::vector<std::array<double, 5>> dumpState();
    // dump the current state of the particles
    // a 2D array of size nParticles x [mass, radius, x, y, speed]
    std::array<double, 4> statistics();
    // return [mean, std] for the current velocity distribution
};
