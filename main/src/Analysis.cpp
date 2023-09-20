#include "main.h"

class Analysis{
  public:
  int particleNum, strands, i;
  double safeMultiplier = 1.4; // Multiplier with safe distance
  std::string type, output,topology;
  LR_vector box, energy;
  std::vector<Particle> particles;
  std::vector<std::vector<Traj>> traj; // For storing trajectory;



  Traj trajtemp;
};