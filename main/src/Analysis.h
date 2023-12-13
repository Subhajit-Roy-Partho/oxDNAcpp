#include "main.h"
#include <random>
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
  int particleNum=0, strands=0, i,tempInt,particlePerStrand=0,particleTypes=1;
  double safeMultiplier = 1.4; // Multiplier with safe distance
  std::string type, output,topology,tempString,config;
  LR_vector box, energy;
  std::vector<Particle> particles;
  std::vector<std::vector<Traj>> traj; // For storing trajectory;

  Traj trajtemp;

  Analysis(std::string topology="", std::string config="", std::string type = "", std::string output = "output", std::string externalForces = "", std::string parameter1 = "", std::string parameter2 = "");
  ~Analysis();
  // Output the center for a number of index
  LR_vector CenterForIndex(int *indexes, int N=-1); // if int* is passed 
  LR_vector CenterForIndex(vector<int> indexes); // Taking advantage of modern c++ if vector is passed.
  LR_vector CenterForIndex(int N);// to get the mean position.
  // template <typename A, typename B> LR_vector center(A indexes, B N);
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
  bool writeCrystalTopology(std::string topology="");
  bool writeDNAtopology(std::string topology="");
  bool writeConfig(std::string config = "");
  bool writeCCGtopology(string topology = "");
  bool writeCCGviewTopology(string topology="");
  bool shiftbox(LR_vector shift={0,0,0}); // Shift the box by given amount but for default case resizes the box
  bool testBoxOverloaded();
  bool generatePSP(Analysis *PSP,vector<vector<int>>ids,vector<int> &colors,vector<double> &radius,int numNeighbour=5,double fixedSpringConstant=0); // Generate PSP particles from DNA particles
  bool addNewType(LR_vector shift ,vector<int> colors={},vector<double> radius={});
  bool addFalseType(vector<int> colors);
  bool populateSingle(int num=125,double seperator=1);
  bool populate(int num=125,double seperator=1);
  bool boxToCubic();
  bool readCCGtopology(std::string topology);
  bool writeMGL(std::string output);
  string line, temp;
  istringstream ss;

  double subBoxing(double coordinate, double divisor);
  bool readCrystalTopology(std::string topology);
  bool readCrystalPatches(std::string patch);
  bool readCrystalParticles(std::string particles);
  bool readDNAtopology(std::string topology);
  bool readConfig(std::string config="");
  bool readPatches(std::string patches);
  bool readParticles(std::string crystalpar);


  bool selectIDs(Analysis* selected,std::vector<int> ids,bool reboxing=true); // select only few particles
  bool reboxing(int offset=2); //optimize the box dimensions

  vector<int> draw();

private:
};


template <typename T> vector<size_t> sort_indexes(const vector<T> &v);