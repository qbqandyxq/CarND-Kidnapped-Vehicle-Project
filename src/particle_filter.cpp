/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

static default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    // Number of particles to draw
<<<<<<< HEAD
    num_particles = 10;
    default_random_engine gen;
=======
    num_particles = 40;

>>>>>>> 26c0254b7e951a8508cb50536bba0c4b81b39e94
    //Standard deviations for x,y, theta
    //create a normal Gaussian noise distribution for x,y, theta
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    //
    for (int i=0; i<num_particles; i++) {

        Particle init_particle;
        init_particle.id = i;
        init_particle.x=dist_x(gen);
        init_particle.y=dist_y(gen);
        init_particle.theta=dist_theta(gen);
        init_particle.weight=1.0;
        particles.push_back(init_particle);
        weights.push_back(1.0);
    }
    is_initialized=true;
}

<<<<<<< HEAD
//void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
//    // TODO: Add measurements to each particle and add random Gaussian noise.
//    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
//    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
//    //  http://www.cplusplus.com/reference/random/default_random_engine/
//    double dx=0.0;
//    double dy=0.0;
//    double dtheta=0.0;
//    default_random_engine gen;
//    for (int i = 0; i<num_particles; i++){
////        Particle &currentParticle = particles[i];
//        if(yaw_rate>0.0001){
//            dx=velocity/yaw_rate*(sin(particles[i].theta + delta_t* yaw_rate) - sin(particles[i].theta));
//            dy=velocity/yaw_rate*(cos(particles[i].theta)- cos(particles[i].theta + yaw_rate*delta_t));
//            dtheta = yaw_rate * delta_t;
//        }else{
//            dx=velocity * cos(particles[i].theta)*delta_t;
//            dy=velocity * sin(particles[i].theta)*delta_t;
////            dtheta = 0.0;
//        }
//        normal_distribution<double> dist_dx(dx, std_pos[0]);
//        normal_distribution<double> dist_dy(dy, std_pos[1]);
//        normal_distribution<double> dist_dtheta(dtheta, std_pos[2]);
//
//        particles[i].x += dist_dx(gen);
//        particles[i].y += dist_dy(gen);
//        particles[i].theta += dist_dtheta(gen);
//    }
//
//}
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
    double dx = 0;
    double dy = 0;
    double dtheta = 0;
    default_random_engine gen;
    for (unsigned int i = 0; i != num_particles; ++i)
    {
        Particle &currentParticle = particles[i];
        
        if (fabs(yaw_rate) > 0.0001)
        {
            dx = velocity/yaw_rate*(sin(currentParticle.theta+yaw_rate*delta_t)-sin(currentParticle.theta));
            dy = velocity/yaw_rate*(cos(currentParticle.theta)-cos(currentParticle.theta+yaw_rate*delta_t));
            dtheta = yaw_rate*delta_t;
        }
        else
        {
            dx = velocity*cos(currentParticle.theta)*delta_t;
            dy = velocity*sin(currentParticle.theta)*delta_t;
            dtheta = 0;
=======
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    double dx=0.0;
    double dy=0.0;
    double dtheta=0.0;
    
    for (int i = 0; i<num_particles; i++){
//        Particle &currentParticle = particles[i];
        if(fabs(yaw_rate)>0.0001){
            dx=velocity/yaw_rate*(sin(particles[i].theta + delta_t* yaw_rate) - sin(particles[i].theta));
            dy=velocity/yaw_rate*(cos(particles[i].theta)- cos(particles[i].theta + yaw_rate*delta_t));
            dtheta = yaw_rate * delta_t;
        }else{
            dx=velocity * cos(particles[i].theta)*delta_t;
            dy=velocity * sin(particles[i].theta)*delta_t;
//            dtheta = 0.0;
>>>>>>> 26c0254b7e951a8508cb50536bba0c4b81b39e94
        }
        normal_distribution<double> dist_dx(dx, std_pos[0]);
        normal_distribution<double> dist_dy(dy, std_pos[1]);
        normal_distribution<double> dist_dtheta(dtheta, std_pos[2]);

        particles[i].x += dist_dx(gen);
        particles[i].y += dist_dy(gen);
        particles[i].theta += dist_dtheta(gen);
    }
}
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    // Nearest neighbor
    for(unsigned int i=0;i<observations.size();i++){
        double mindistance = 99999;
        int match_index = 0;
        for(unsigned int j=0;j<predicted.size();j++){
            double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            if(distance < mindistance){
                mindistance=distance;
                match_index = j;
            }
        }
        observations[i].id = match_index;
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    for (unsigned int i = 0; i < num_particles; ++i)
    {
        Particle &currentParticle = particles[i];
        std::vector<LandmarkObs> predicted;
        std::vector<LandmarkObs> observationsT = observations;
        
        for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j)
        {
            LandmarkObs predictedLandMark;
            predictedLandMark.x = map_landmarks.landmark_list[j].x_f;
            predictedLandMark.y = map_landmarks.landmark_list[j].y_f;
            predictedLandMark.id = map_landmarks.landmark_list[j].id_i;
            
            predicted.push_back(predictedLandMark);
        }
        
        for (unsigned int j = 0; j< observationsT.size(); ++j)
        {
            
            observationsT[j].x = currentParticle.x + (cos(currentParticle.theta)*observations[j].x) - (sin(currentParticle.theta) * observations[j].y);
            observationsT[j].y = currentParticle.y + (sin(currentParticle.theta)*observations[j].x) + (cos(currentParticle.theta) * observations[j].y);
        }
        
        dataAssociation(predicted, observationsT);
        
        double weight = 1.0;
        for (unsigned int j = 0; j< observationsT.size(); ++j)
        {
            double sig_x = std_landmark[0];
            double sig_y = std_landmark[1];
            double x_obs = observationsT[j].x;
            double y_obs = observationsT[j].y;
            double mu_x = map_landmarks.landmark_list[observationsT[j].id].x_f;
            double mu_y = map_landmarks.landmark_list[observationsT[j].id].y_f;
            
            // calculate normalization term
            double gauss_norm = (1/(2 * M_PI * sig_x * sig_y));
            
            // calculate exponent
            double exponent = pow(x_obs - mu_x, 2)/(2*pow(sig_x, 2))
            + pow(y_obs - mu_y, 2)/(2*pow(sig_y, 2));
            
            // calculate weight using normalization terms and exponent
            weight *= gauss_norm * exp(-exponent);
        }
        
        currentParticle.weight = weight;
        weights[i] = weight;
    }
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    

    std::vector<Particle> newParticles;
    std::discrete_distribution<int> weightsDistribution (weights.begin(), weights.end());
    for (unsigned int i = 0; i < num_particles; ++i)
    {
        newParticles.push_back(particles[weightsDistribution(gen)]);
    }
    particles = newParticles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    
    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();
    
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
