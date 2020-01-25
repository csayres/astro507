#include <iostream>
#include <cmath>
#include "particles.h"

double randomSample(){
    // return between 0 and 1
    return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}


Particle::Particle(double mass, double radius, double x, double y, double vx, double vy)
    : mass(mass), radius(radius), x(x), y(y), vx(vx), vy(vy){};


std::array<double, 2> Particle::nextPosition(double dt){
    std::array<double, 2> nextXY;
    nextXY[0] = x + vx*dt;
    nextXY[1] = y + vy*dt;
    return nextXY;

}

void Particle::updatePosition(double dt){
    auto nextXY = nextPosition(dt);
    x = nextXY[0];
    y = nextXY[1];
}

void Particle::setVelocities(double vxnext, double vynext){
    vx = vxnext;
    vy = vynext;
}

double Particle::speed(){
    return hypot(vx, vy);
}

double Particle::direction(){
    return atan2(vy, vx);
}


Box::Box(double width, double height, int seed)
    : width(width), height(height)
{
    srand(seed);
}

std::array<double, 2> Box::randomPoint(double radius){
    std::array<double, 2> xy;
    xy[0] = randomSample()*(width - 2*radius) + radius;
    xy[1] = randomSample()*(height - 2*radius) + radius;
    return xy;
}

bool Box::isCollided(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
    double dx, dy, dist, nextDist;
    std::array<double, 2> p1next, p2next;

    dx = p1->x - p2->x;
    dy = p1->y - p2->y;
    dist = hypot(dx, dy); // current distance

    // first check if particles are moving towards eachother
    // if not, don't consider them collided
    p1next = p1->nextPosition(dt);
    p2next = p2->nextPosition(dt);
    nextDist = hypot(p1next[0]-p2next[0], p1next[1]-p2next[1]);
    if (nextDist > dist){
        // particles are moving away (collision probably already happened)
        return false;
    }

    if (dist <= p1->radius + p2->radius){
        // distance is within the collision zone
        return true;
    }

    return false;
}

void Box::addRandomParticle(double mass, double radius, double speed){
    bool collided;
    double ang, vx, vy;

    ang = randomSample()*2*M_PI;
    vx = cos(ang)*speed;
    vy = sin(ang)*speed;

    if (particles.empty()){
        // this is the first particle added
        // add it and exit
        auto xy = randomPoint(radius);
        auto p = std::make_shared<Particle>(mass, radius, xy[0], xy[1], vx, vy);
        particles.push_back(p);
        return;
    }

    // ensure we find a non-colliding position for the particle
    // try 100000 times before giving up and throwing a runtime error
    int maxIter = 100000;
    for (int ii=0; ii<=maxIter; ii++){
        collided = false;
        auto xy = randomPoint(radius);
        auto p = std::make_shared<Particle>(mass, radius, xy[0], xy[1], vx, vy);
        for (auto p2 : particles){
            if (isCollided(p, p2)){
                collided = true;
                break;
            }
        }

        if (!collided){
            particles.push_back(p);
            std::cout << "particles " << particles.size() << std::endl;
            break;
        }
        if (ii==maxIter){
            throw std::runtime_error("cannot find place for new particle");
        }
    }
}

void Box::handleCollision(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2){
    // assumes collision has been verified for input particle indices...
    // implement the equations from:
    // https://en.wikipedia.org/wiki/Elastic_collision
    double contactAngle, A1, A2, v1x, v1y, v2x, v2y;

    contactAngle = atan2(
        p1->y - p2->y,
        p1->x - p2->x
    );

    A1 = (
        p1->speed() *
        cos(p1->direction() - contactAngle) *
        (p1->mass - p2->mass) +
        2 * p2->mass * p2->speed() *
        cos(p2->direction() - contactAngle)
    ) / (p1->mass + p2->mass);

    A2 = (
        p2->speed() *
        cos(p2->direction() - contactAngle) *
        (p2->mass - p1->mass) +
        2 * p1->mass * p1->speed() *
        cos(p1->direction() - contactAngle)
    ) / (p2->mass + p1->mass);

    // new x velocity for particle 1
    v1x = A1 * cos(contactAngle) +
        p1->speed() * sin(p1->direction() - contactAngle) *
        cos(contactAngle + M_PI / 2);
    // new y velocity for particle 1
    v1y = A1 * sin(contactAngle) + \
        p1->speed() * sin(p1->direction() - contactAngle) *
        sin(contactAngle + M_PI / 2);

    // new x velocity for particle 2
    v2x = A2 * cos(contactAngle) + \
        p2->speed() * sin(p2->direction() - contactAngle) *
        cos(contactAngle + M_PI / 2);
    // new y velocity for particle 2
    v2y = A2 * sin(contactAngle) + \
        p2->speed() * sin(p2->direction() - contactAngle) *
        sin(contactAngle + M_PI / 2);

    p1->setVelocities(v1x, v1y);
    p2->setVelocities(v2x, v2y);

}

std::vector<std::array<double, 5>> Box::dumpState(){
    std::vector<std::array<double, 5>> ensemble;

    for (auto p : particles){
        std::array<double, 5> pState;
        pState[0] = p->mass;
        pState[1] = p->radius;
        pState[2] = p->x;
        pState[3] = p->y;
        pState[4] = p->speed();
        ensemble.push_back(pState);
    }
    return ensemble;
}

void Box::runSim(int steps, double timestep, int saveEvery){
    double pressure, collisions;

    dt = timestep;
    int nParticles = particles.size();

    for (int step=0; step < steps; step++){
        // begin main sim loop
        // first propatate particles
        pressure = 0;
        collisions = 0;

        for (auto p : particles){
            p->updatePosition(dt);
            // reflect the particle if it's at the boundary
            if ( (p->x + p->radius >= width) and (p->vx > 0)){
                pressure += p->mass * 2 * p->vx / height;
                p->setVelocities(-1*p->vx, p->vy);
            }
            if ( (p->y + p->radius >= height) and (p->vy > 0)){
                p->setVelocities(p->vx, -1*p->vy);
                pressure += p->mass * 2 * p->vy / width;
            }
            if ( (p->x - p->radius <= 0) and (p->vx < 0)){
                p->setVelocities(-1*p->vx, p->vy);
                pressure += p->mass * 2 * p->vx / height;
            }
            if ( (p->y - p->radius <= 0) and (p->vy < 0)){
                p->setVelocities(p->vx, -1*p->vy);
                pressure += p->mass * 2 * p->vy / width;
            }
        }
        pressureSteps.push_back(pressure);


        // next check all pairwise combos of particles
        // for collisions, update velocities if necessary
        // upper triangular nested loop
        for (int ii=0; ii < nParticles - 1; ii++){
            for (int jj=ii+1; jj < nParticles; jj++){
                auto p1 = particles[ii];
                auto p2 = particles[jj];
                if (isCollided(p1, p2)){
                    handleCollision(p1, p2);
                    collisions += 1;
                }
            }
        }
        collisionSteps.push_back(collisions);


        if ((step % saveEvery) == 0){
            double fraction = float(step)/(float)steps;
            std::cout << "saving at step " << step << " fraction done " << fraction << std::endl;
            particleSteps.push_back(dumpState());
        }

    }

}







