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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).	
	// generator	
	default_random_engine gen;	
	num_particles = 200;

	weights.resize(num_particles);
	particles.resize(num_particles);
	double std_x, std_y, std_theta, sample_x,sample_y,sample_theta;
	std_x = 2.0f;
	std_y = 2.0f;
	std_theta = 0.05f;

	// generating gaussians
	normal_distribution<double> dist_x(x,std_x);
	normal_distribution<double> dist_y(y,std_y);
	normal_distribution<double> dist_theta(theta,std_theta);

	for ( int i = 0; i<num_particles; i++) {
		particles[i].id=i;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
		// inintial weight 1
		particles[i].weight = 1.0;	
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//performing calculations for predictions
	const float theta_prime = yaw_rate * delta_t;
	
	for (int i=0; i<num_particles; i++) {
	const float theta = particles[i].theta;
	
	if(abs(yaw_rate)>.00001){
		particles[i].x += (velocity/yaw_rate) * (sin(theta + theta_prime) - sin(theta));
		particles[i].y += (velocity/yaw_rate) * (cos(theta) - cos(theta + theta_prime));
		particles[i].theta += yaw_rate * delta_t;
	} // end if

		// calculation different with nonzero-yaw_rate avoid divide by zero
		// also yaw_rate being zero means no change to theta 
	else {
		particles[i].x += velocity * delta_t * cos(theta);
		particles[i].y += velocity * delta_t * sin(theta);	
	} // end else
	
	// adding noise

	// generator
	default_random_engine gen;
	const float std_x = 2.0f;
	const float std_y = 2.0f;
	const float std_theta = 0.05f;

	// generating gaussians for random noise, mean zero
	normal_distribution<double> dist_x(0,std_x);
	normal_distribution<double> dist_y(0,std_y);
	normal_distribution<double> dist_theta(0,std_theta);

	particles[i].x += dist_x(gen);
	particles[i].y += dist_y(gen);
	particles[i].theta += dist_theta(gen);
	
	} // end forloop i
	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// determine size length of landmarks
	/*	
	int num_observations = observations.size();
	
	for (int i=0; i<observations.size(); i++) {
		// particle  to map coordinate transformations:
		obj_x = observations[i].x	
		
	} // end forloop i

	*/

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
	// std_landmark Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]

	const double constant_term = 1/ ( 2* M_PI * std_landmark[0] * std_landmark[1]);
	const double sigma_x = 2 * std_landmark[0] * std_landmark[0];
	const double sigma_y = 2 * std_landmark[1] * std_landmark[1];

	// loop through particles
	for (int i=0; i<num_particles; i++) {
	// initializing multivariate gaussian distribution for particles
	double gaussian = 1.0f;
	
	// caching values for individual particles 
	double sin_term = sin(particles[i].theta);
	double cos_term = cos(particles[i].theta);
	
		// looping through each observation of a particle[i]	
		for (int j=0; j < observations.size(); j++) {

		double x_trans,y_trans ;
		x_trans = observations[j].x * cos_term - observations[j].y * sin_term + particles[i].x;
		y_trans = observations[j].x * sin_term + observations[j].y * cos_term + particles[i].y;

		// nearest neighbor for closest landmark 
		vector<Map:: single_landmark_s> landmarks = map_landmarks.landmark_list;
		vector<double> landmark_dist (landmarks.size()); // size optimization

			// looping through landmarks list to find closest one to each transformed observation
			for(int k = 0 ; k < landmarks.size() ; k++ ) {
			// using sensor_range to look for proximal landmarks only
			const double landmark_indiv_dist = dist(particles[i].x,particles[i].y,landmarks[k].x_f,landmarks[k].y_f );
	//sqrt(pow(particles[i].x - landmarks[k].x_f,2) + pow(particles[i].y - landmarks[k].y_f, 2));
				
				if(landmark_indiv_dist <= sensor_range) {
				landmark_dist[k] = sqrt( pow(x_trans - landmarks[k].x_f,2) + pow(y_trans-landmarks[k].y_f,2));
				} // end if
				else {
				// giving maximum distance to landmarks far from sensors, excluding from consideration
				landmark_dist[k] = 1000000.0;

				} // end else
			} // end forloop k
		// landmark to observation association 
		int min_index = distance(landmark_dist.begin(),min_element(landmark_dist.begin(),landmark_dist.end()));
		
		double nx = landmarks[min_index].x_f;	//landmarks associated within sensor range
		double ny = landmarks[min_index].y_f;
		// calculating dx,dy
		double dx = x_trans - nx;
		double dy = y_trans - ny;
		double exp_term = ((dx*dx)/sigma_x) + ((dy*dy)/sigma_y);
		
		// multiplying to gaussian		
		gaussian *= constant_term * exp(-exp_term);

		} // end forloop j
	// add particle weights for each particle
	particles[i].weight = gaussian;
	weights[i] = particles[i].weight; // for total weight normalization

	} // end forloop i 
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// new vector of particles
	vector<Particle> resampled_particles (num_particles);
	default_random_engine gen;
	// looping through particles
	for (int i=0; i<num_particles; i++) {
		std::discrete_distribution<int> index(weights.begin(),weights.end());
		resampled_particles[i] = particles[index(gen)];	

	} // end forloop i
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
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
