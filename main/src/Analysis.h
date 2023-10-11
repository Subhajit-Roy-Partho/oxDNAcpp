#include "main.h"
#include <eigen3/Eigen/Dense>
using namespace std;

class Forces
{
public:
  std::vector<int> particles;
  std::vector<float> stiff, rate, position;
  std::vector<std::string> type;
  std::vector<LR_vector> dir, pos0;

  bool add(string name, int particle = -1, float stiff = 0, LR_vector normal = {0, 0, 1}, LR_vector pos = {0, 0, 0});
  bool addRepulsion(int particles = -1, float stiff = 1, LR_vector dir = {0, 0, 1}, float position = 0);
  bool addHarmonic(int particles = 0, float stiff = 1, float rate = 0, LR_vector pos0 = {0, 0, 0}, LR_vector dir = {1, 0, 0});
  bool addHarmonicFromParticle(Particle *particle, float stiff = 1, float rate = 0);
  bool save(string filename = "external_forces.txt");
};


class Analysis{
  public:
  int particleNum, strands, i,tempInt;
  double safeMultiplier = 1.4; // Multiplier with safe distance
  std::string type, output,topology,tempString;
  LR_vector box, energy;
  std::vector<Particle> particles;
  std::vector<std::vector<Traj>> traj; // For storing trajectory;

  Traj trajtemp;

  Analysis(std::string topology, std::string config, std::string type = "", std::string output = "output", std::string externalForces = "", std::string parameter1 = "", std::string parameter2 = "");
  ~Analysis();
  // Output the center for a number of index
  LR_vector CenterForIndex(int *indexes, int N);
  LR_vector NormalToPlane(Eigen::MatrixXd points);
  void inboxing(LR_vector center = {0, 0, 0});
  bool customSeedForces(std::vector<int> ids);
  bool correctA(string newFile);
  bool randomReplaceColor(int originalColor, int newColor, vector<int> ignore = {}, int N = 1);
  bool randomReplacePosition(int originalColor, Particle *newParticle, vector<int> ignore = {});
  bool pickAndPlace(int *cluster, int N, string tname, LR_vector centralShift = {0, 0, 0});

  // Computes the minimum distances for two particles with indices p and q
  template <typename A, typename B>double min_image(A p, B q);
  bool planeFitting(int *cluster, LR_vector center, LR_vector normal);
  bool writeCrystalTopology(std::string topology = "");
  bool writeDNAtopology(std::string topology);
  bool writeConfig(std::string config = "");
  bool shiftbox(LR_vector shift={0,0,0});
  bool testBoxOverloaded();
  bool generatePSP(Analysis *PSP,vector<vector<int>>,int numCluster=12,int avgSize=10,int numNeighbour=5); // Generate PSP particles from DNA particles

private:
  string line, temp;
  istringstream ss;

  double subBoxing(double coordinate, double divisor);
  int readCrystalTopology(std::string topology);
  bool readDNAtopology(std::string topology);
  int readConfig(std::string config);
  int readPatches(std::string patches);
  int readParticles(std::string crystalpar);
};


template <typename T> vector<size_t> sort_indexes(const vector<T> &v);