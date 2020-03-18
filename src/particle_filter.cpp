/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::discrete_distribution;
using std::default_random_engine;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  unsigned int num_particles = 100;  // TODO: Set the number of particles
  
  default_random_engine gen;
  normal_distribution<> dist_x{x, std[0]};
  normal_distribution<> dist_y{y, std[1]};
  normal_distribution<> dist_theta{theta, std[2]};
  
  for (unsigned int i = 0; i < num_particles; ++i) {
  	Particle p;
    
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1;
    
    particles.push_back(p);
  }

  is_initialized = true;

  std::cout << "The Particle Filter is initialized with " << num_particles << " particles!" << std::endl;
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
  default_random_engine gen;
  normal_distribution<> dist_x_noise{0, std_pos[0]};
  normal_distribution<> dist_y_noise{0, std_pos[1]};
  normal_distribution<> dist_theta_noise{0, std_pos[2]};
  
  if (yaw_rate > 0.0001 || yaw_rate < -0.0001) {
    for (Particle &p : particles) {
      p.x += (velocity / yaw_rate) * (sin(p.theta + yaw_rate * delta_t) - sin(p.theta)) + dist_x_noise(gen);
      p.y += (velocity / yaw_rate) * (cos(p.theta) - cos(p.theta + yaw_rate * delta_t)) + dist_y_noise(gen);
      p.theta += yaw_rate * delta_t + dist_theta_noise(gen);
    }
  } else {
    for (Particle &p : particles) {
      p.x += velocity * delta_t * cos(p.theta) + dist_x_noise(gen);
      p.y += velocity * delta_t * sin(p.theta) + dist_y_noise(gen);
    }
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
  weights.clear();
  for (Particle &p : particles) {
    // reset particles variables
    p.weight = 1;
    p.associations.clear();
    p.sense_x.clear();
    p.sense_y.clear();
    
    // select only the landmarks in the sensor range
    vector<Map::single_landmark_s> landmarks_in_range;
    findLandmarksInRange(landmarks_in_range, map_landmarks, p, sensor_range);
    
    for (auto observation: observations) {
      // Convert observation coordinates to map coordinates
      LandmarkObs map_observation;
      convert2MapCoordinates(map_observation, observation, p);
      
      // associate a landmark to each observation and update the particle variables
      Map::single_landmark_s associated_landmark = dataAssociation(landmarks_in_range, map_observation);
      p.associations.push_back(map_observation.id);
      p.sense_x.push_back(map_observation.x);
      p.sense_y.push_back(map_observation.y);
      
      // calculate the weight of the particle
      double observation_probability = 
        1 / (2 * M_PI * std_landmark[0] * std_landmark[1]) * 
        exp(-(
               (pow(map_observation.x - associated_landmark.x_f, 2) / (2 * pow(std_landmark[0], 2))) + 
               (pow(map_observation.y - associated_landmark.y_f, 2) / (2 * pow(std_landmark[1], 2)))
            ) );
      p.weight *= observation_probability; 
    }
    weights.push_back(p.weight);
  }
}

/**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
void ParticleFilter::resample() {
  default_random_engine gen;
  discrete_distribution<int> ddist(weights.begin(), weights.end());
  
  vector<Particle> new_particles;
  for (unsigned int i = 0; i < particles.size(); ++i) {
    new_particles.push_back(particles[ddist(gen)]);
  }
  particles = new_particles;
}

// Helper functions
void ParticleFilter::findLandmarksInRange(std::vector<Map::single_landmark_s> &landmarks_in_range, const Map &map_landmarks, const Particle &p, double sensor_range) {
  for (auto lm : map_landmarks.landmark_list) {
    double distance = dist(p.x, p.y, lm.x_f, lm.y_f);
    
    if (distance <= sensor_range) {
      Map::single_landmark_s landmark_in_range;
      landmark_in_range.id_i = lm.id_i;
      landmark_in_range.x_f = lm.x_f;
      landmark_in_range.y_f = lm.y_f;      
      landmarks_in_range.push_back(landmark_in_range);
      }
    }  
}

void ParticleFilter::convert2MapCoordinates(LandmarkObs &map_observation, const LandmarkObs &observation, const Particle &p) {
  map_observation.x = p.x + cos(p.theta) * observation.x - sin(p.theta) * observation.y;
  map_observation.y = p.y + sin(p.theta) * observation.x + cos(p.theta) * observation.y;
}

Map::single_landmark_s ParticleFilter::dataAssociation(vector<Map::single_landmark_s> landmarks, 
                                     LandmarkObs &observation) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  double min_dist = dist(observation.x, observation.y, landmarks[0].x_f, landmarks[0].y_f);
  Map::single_landmark_s associated_landmark = landmarks[0];
  
  observation.id = landmarks[0].id_i;
  for(auto landmark : landmarks) {
    double temp_dist = dist(observation.x, observation.y, landmark.x_f, landmark.y_f);

    if (temp_dist < min_dist) {
      min_dist = temp_dist;
      associated_landmark = landmark;
      observation.id = landmark.id_i;
    }
  }
  
  return associated_landmark;
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