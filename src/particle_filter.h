/**
 * particle_filter.h
 * 2D particle filter class.
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <string>
#include <vector>
#include "helper_functions.h"

struct Particle {
  int id;
  double x;
  double y;
  double theta;
  double weight;
  std::vector<int> associations;
  std::vector<double> sense_x;
  std::vector<double> sense_y;
};


class ParticleFilter {  
 public:
  // Constructor
  // @param num_particles Number of particles
  ParticleFilter() : is_initialized(false) {}

  // Destructor
  ~ParticleFilter() {}

  /**
   * init Initializes particle filter by initializing particles to Gaussian
   *   distribution around first position and all the weights to 1.
   * @param x Initial x position [m] (simulated estimate from GPS)
   * @param y Initial y position [m]
   * @param theta Initial orientation [rad]
   * @param std[] Array of dimension 3 [standard deviation of x [m], 
   *   standard deviation of y [m], standard deviation of yaw [rad]]
   */
  void init(double x, double y, double theta, double std[]);

  /**
   * prediction Predicts the state for the next time step
   *   using the process model.
   * @param delta_t Time between time step t and t+1 in measurements [s]
   * @param std_pos[] Array of dimension 3 [standard deviation of x [m], 
   *   standard deviation of y [m], standard deviation of yaw [rad]]
   * @param velocity Velocity of car from t to t+1 [m/s]
   * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
   */
  void prediction(double delta_t, double std_pos[], double velocity, 
                  double yaw_rate);
  
  /**
   * updateWeights Updates the weights for each particle based on the likelihood
   *   of the observed measurements. 
   * @param sensor_range Range [m] of sensor
   * @param std_landmark[] Array of dimension 2
   *   [Landmark measurement uncertainty [x [m], y [m]]]
   * @param observations Vector of landmark observations
   * @param map Map class containing map landmarks
   */
  void updateWeights(double sensor_range, double std_landmark[], 
                     const std::vector<LandmarkObs> &observations,
                     const Map &map_landmarks);
  
  /**
   * resample Resamples from the updated set of particles to form
   *   the new set of particles.
   */
  void resample();

  /**
   * Used for obtaining debugging information related to particles.
   */
  std::string getAssociations(Particle best);
  std::string getSenseCoord(Particle best, std::string coord);
  
  /**
   * initialized Returns whether particle filter is initialized yet or not.
   */
  bool initialized() const {
    return is_initialized;
  }

  // Set of current particles
  std::vector<Particle> particles;

 private:
  // Flag, if filter is initialized
  bool is_initialized;
  
  // Vector of weights of all particles
  std::vector<double> weights;
  
  // Find the landmarks in the particle range
  void findLandmarksInRange(std::vector<Map::single_landmark_s> &landmarks_in_range, 
                            const Map &map_landmarks, const Particle &p, double sensor_range);
  
  // Convert local coordinates to map coordinates
  void convert2MapCoordinates(LandmarkObs &map_observation, const LandmarkObs &observation, 
                              const Particle &p);
  
  /**
   * dataAssociation Finds which landmark corresponds to a provided observation 
   *   (likely by using a nearest-neighbors data association).
   * @param landmarks vector
   * @param observation
   * @return associated landmark
   */
  Map::single_landmark_s dataAssociation(std::vector<Map::single_landmark_s> landmarks, 
                         LandmarkObs& observation);
};

#endif  // PARTICLE_FILTER_H_