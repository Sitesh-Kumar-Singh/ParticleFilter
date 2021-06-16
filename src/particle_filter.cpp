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
#include "helper_functions.h"

using std::string;
using std::vector;
;
static std::default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  
  // init particles
  for (int i = 0; i < num_particles; i++) {
    Particle p_init;
    p_init.id = i;
    p_init.x = dist_x(gen);
    p_init.y = dist_y(gen);
    p_init.theta = dist_theta(gen);
    p_init.weight = 1.0;

    particles.push_back(p_init);
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // define normal distributions for sensor noise
   std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);
  
  for (int i = 0; i < num_particles; i++) {
    // calculate new state
    if (fabs(yaw_rate) < 0.00001) {  
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    } 
    else {
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    // add noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_theta(gen);
    particles[i].theta += dist_theta(gen);
  }


}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (int i=0;i<observations.size();i++){
    LandmarkObs obs = observations[i];
    double min_dist = std::numeric_limits<double>::max();
    int map_id = -1;
    for (int j=0; j<predicted.size();j++){
      LandmarkObs current_pre = predicted[j];
      double distance = dist(obs.x,obs.y,current_pre.x,current_pre.y);
      if (distance < min_dist){
        min_dist = distance;
        map_id = current_pre.id;
      }
    }
    observations[i].id = map_id;
  }

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   /*
    1. get x,y of particles.
     2. get x,y of landmark in sensor range
     3. transfrom vehicle co-ordinate to map co-ordinates
     4. perform data association
     5. Calculate weight with multivariate gaussian
     */
  // Step 1
   for (int i = 0; i < num_particles; i++) {

        double p_x = particles[i].x;
        double p_y = particles[i].y;
        double p_theta = particles[i].theta;
     
     
     //Step 2
        vector<LandmarkObs> predictions;
        for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

            float landmark_x = map_landmarks.landmark_list[j].x_f;
            float landmark_y = map_landmarks.landmark_list[j].y_f;
            int landmark_id = map_landmarks.landmark_list[j].id_i;

            if (fabs(landmark_x - p_x) <= sensor_range && fabs(landmark_y - p_y) <= sensor_range) {

                 predictions.push_back(LandmarkObs{ landmark_id, landmark_x, landmark_y });
            }
        }
       
       //Step 3
        vector<LandmarkObs> transf_obs; 
        for (unsigned int k = 0; k < observations.size(); k++) {
            double t_x = cos(p_theta)*observations[k].x - sin(p_theta)*observations[k].y + p_x;
            double t_y = sin(p_theta)*observations[k].x + cos(p_theta)*observations[k].y + p_y;
            transf_obs.push_back(LandmarkObs{ observations[k].id, t_x, t_y });
        }
    
        dataAssociation(predictions, transf_obs);

        particles[i].weight = 1.0;
        for (unsigned int l = 0; l < transf_obs.size(); l++) {

            double obs_x, obs_y;
            double pre_x, pre_y;
            obs_x = transf_obs[l].x;
            obs_y = transf_obs[l].y;
            int associated_prediction = transf_obs[l].id;

            for (unsigned int m = 0; m < predictions.size(); m++) {
                if (predictions[m].id == associated_prediction) {
                    pre_x = predictions[m].x;
                    pre_y = predictions[m].y;
                }
            }
            double std_x = std_landmark[0];
            double std_y = std_landmark[1];
            double obs_weight = ( 1/(2*M_PI*std_x*std_y)) * exp( -( pow(pre_x-obs_x,2)/(2*pow(std_x, 2)) + (pow(pre_y-obs_y,2)/(2*pow(std_y, 2))) ) );

            particles[i].weight *= obs_weight;
        }
       
    
    }
     
}


void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
vector<Particle> sampled_particles;

  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }
  std::uniform_int_distribution<int> uni_distr_index(0,num_particles -1);
  auto rnd_index = uni_distr_index(gen);
  
  double max_weight = *max_element(weights.begin(),weights.end());
  
  std::uniform_real_distribution<double> uni_dist_weight(0.0, 2*max_weight);
  
  double beta = 0.0;
  
  
  for (int i = 0; i < num_particles; i++) {
    beta += uni_dist_weight(gen);
    while (beta > weights[rnd_index]) {
      beta -= weights[rnd_index];
      rnd_index = (rnd_index + 1) % num_particles;
    }
    sampled_particles.push_back(particles[rnd_index]);
  }
  particles = sampled_particles;
  
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}