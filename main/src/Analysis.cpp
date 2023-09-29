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

  bool add(string name, int particle = -1, float stiff = 0, LR_vector normal = {0, 0, 1}, LR_vector pos = {0, 0, 0})
  {
    type.push_back(name);
    particles.push_back(particle);
    dir.push_back(normal);
    pos0.push_back(pos);
    this->stiff.push_back(stiff);
    return true;
  }
  bool addRepulsion(int particles = -1, float stiff = 1, LR_vector dir = {0, 0, 1}, float position = 0)
  {
    type.push_back("repulsion_plane");
    this->particles.push_back(particles);
    this->stiff.push_back(stiff);
    this->dir.push_back(dir);
    this->position.push_back(position);
    return true;
  }
  bool addHarmonic(int particles = 0, float stiff = 1, float rate = 0, LR_vector pos0 = {0, 0, 0}, LR_vector dir = {1, 0, 0})
  {
    type.push_back("trap");
    this->particles.push_back(particles);
    this->stiff.push_back(stiff);
    this->rate.push_back(rate);
    this->pos0.push_back(pos0);
    this->dir.push_back(dir);
    return true;
  }

  bool addHarmonicFromParticle(Particle *particle, float stiff = 1, float rate = 0)
  {
    return addHarmonic(particle->id, stiff, rate, particle->r);
  }

  bool save(string filename = "external_forces.txt")
  {
    ofstream file(filename);
    if (!file.is_open())
      return false;
    for (int i = 0; i < type.size(); i++)
    {
      file << "{" << std::endl;
      if (type[i] == "trap")
      {
        file << "type = trap\n";
        file << "particle = " << particles[0] << endl;
        particles.erase(particles.begin());
        file << "pos0 = " << pos0[0].x<<","<<pos0[0].y<<","<<pos0[0].z << endl;
        pos0.erase(pos0.begin());
        file << "stiff = " << stiff[0] << endl;
        stiff.erase(stiff.begin());
        file << "rate = " << rate[0] << endl;
        rate.erase(rate.begin());
        file << "dir = " << dir[0].x<<","<<dir[0].y<<","<<dir[0].z << endl;
        dir.erase(dir.begin());
      }
      else if (type[i] == "repulsion_plane")
      {
        file << "type = repulsion_plane" << endl;
        file << "particle = " << particles[0] << endl;
        particles.erase(particles.begin());
        file << "stiff = " << stiff[0] << endl;
        stiff.erase(stiff.begin());
        file << "dir = " << dir[0].x<<","<<dir[0].y<<","<<dir[0].z << endl;
        dir.erase(dir.begin());
        file << "position = " << position[0] << endl;
        position.erase(position.begin());
      }
      file << "}\n";
    }
    file.close();
    return true;
  }
};


class Analysis{
  public:
  int particleNum, strands, i;
  double safeMultiplier = 1.4; // Multiplier with safe distance
  std::string type, output,topology;
  LR_vector box, energy;
  std::vector<Particle> particles;
  std::vector<std::vector<Traj>> traj; // For storing trajectory;

  Traj trajtemp;

  Analysis(std::string topology, std::string config, std::string type = "", std::string output = "output", std::string externalForces = "", std::string parameter1 = "", std::string parameter2 = ""){
    if (type == "crystal")
    {
      this->type = type;
      this->output = output;
      this->topology = topology;
      readCrystalTopology(topology);
      readConfig(config);
      inboxing();
      // pickAndPlace();
    }
  }
  ~Analysis()
  {
    // gsl_matrix_free(V);
    // gsl_vector_free(S);
    // gsl_vector_free(work);
  }
  // Output the center for a number of index
  LR_vector CenterForIndex(int *indexes, int N){
    LR_vector mean = {0, 0, 0};
    for (int i = 0; i < N; i++){
      mean += particles[indexes[i]].r;
    }
    mean /= N;
    return mean;
  }

  LR_vector NormalToPlane(Eigen::MatrixXd points)
  {
    Eigen::Matrix3d m(3, 3);
    m << points.transpose() * points;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
    m = svd.matrixV();
    return (LR_vector){m(0, 2), m(1, 2), m(2, 1)};
  };
  void inboxing(LR_vector center = {0, 0, 0})
  {
    for (int i = 0; i < particleNum; i++)
    {
      particles[i].r.x = subBoxing(particles[i].r.x, box.x);
      particles[i].r.y = subBoxing(particles[i].r.y, box.y);
      particles[i].r.z = subBoxing(particles[i].r.z, box.z);
    }
  }

  bool customSeedForces(std::vector<int> ids)
  {
    Forces force;
    force.addRepulsion();
    for (int i = 0; i < ids.size(); i++)
    {
      force.addHarmonicFromParticle(&particles[ids[i]]);
    }
    force.save();
    return true;
  };

  bool correctA(string newFile){
    Analysis newPar(topology,newFile,"crystal");
    for(int i=0;i<particleNum;i++){
      particles[i].a1=newPar.particles[i].a1;
      particles[i].a3=newPar.particles[i].a3;
    }
    return true;
  }

  bool randomReplaceColor(int originalColor, int newColor, vector<int> ignore = {}, int N = 1)
  {
    for (int i = 0; i < N; i++)
    {
      bool found = false;
      int random = rand() % particleNum;
      for (int j = random; j < particleNum; j++)
      {
        if (ignore.size() > 0 && std::find(ignore.begin(), ignore.end(), particles[j].id) != ignore.cend())
          continue;
        if (particles[j].strand == originalColor)
        {
          particles[j].strand = newColor;
          return true;
        }
      }
      if (!found)
      {
        for (int j = random; j >= 0; j--)
        {
          if (ignore.size() > 0 && std::find(ignore.begin(), ignore.end(), particles[j].id) != ignore.cend())
            continue;
          if (particles[j].strand == originalColor)
          {
            particles[j].strand = newColor;
            return true;
          }
        }
      }
      if (!found)
        return false;
    }
    return true;
  }

  bool randomReplacePosition(int originalColor, Particle *newParticle, vector<int> ignore = {})
  {
    srand(time(0));
    bool found = false;
    int random = rand() % particleNum;
    for (int j = random; j < particleNum; j++)
    {
      if (particles[j].strand == originalColor)
      {
        particles[j].r = newParticle->r;
        particles[j].a1 = newParticle->a1;
        particles[j].a3 = newParticle->a3;
        return true;
      }
    }
    for (int j = random; j >= 0; j--)
    {
      if (particles[j].strand == originalColor)
      {
        particles[j].r = newParticle->r;
        particles[j].a1 = newParticle->a1;
        particles[j].a3 = newParticle->a3;
        return true;
      }
    }
    return false;
  }
  bool pickAndPlace(int *cluster, int N, string tname, LR_vector centralShift = {0, 0, 0})
  {
    Analysis target(topology, "BIG.conf", "crystal");
    target.inboxing();
    Eigen::MatrixXd points(N, 3);

    LR_vector center = CenterForIndex(cluster, N);
    double t = std::remainder((double)-7.1, (double)3);
    for (int i = 0; i < N; i++)
    {
      LR_vector some = particles[cluster[i]].r - particles[cluster[i + 1]].r;
      // cout<<some.module()<<endl;
      particles[cluster[i]].r -= center + centralShift;
      points(i, 0) = particles[cluster[i]].r.x;
      points(i, 1) = particles[cluster[i]].r.y;
      points(i, 2) = particles[cluster[i]].r.z;
    }
    inboxing();

    LR_vector normal = NormalToPlane(points);
    // cout<< normal<<endl;

    Eigen::Matrix3d Rot;
    Rot = Eigen::AngleAxisd(acos(normal.x) - M_PI, Eigen::Vector3d::UnitY()) *
          Eigen::AngleAxisd(M_PI - acos(normal.y), Eigen::Vector3d::UnitX());

    // cout<<(LR_vector) {points(0,0),points(0,1),points(0,2)}<<endl;
    points = points * Rot;
    double safeDistance = 99999999999999999;
    double dist;
    for (int i = 0; i < N - 1; i++)
    {
      dist = npNorm((LR_vector){points(i, 0), points(i, 1), points(i, 2)} - (LR_vector){points(i + 1, 0), points(i + 1, 1), points(i + 1, 2)});
      if (dist < safeDistance)
        safeDistance = dist;
    }
    safeDistance *= safeMultiplier;
    // writeConfig();
    // std::cout <<safeDistance<<std::endl;
    std::vector<int> infected, infectedColors;
    // cout <<min_image(1,particles[2].r)<<endl;
    for (int j = 0; j < N; j++)
    {
      for (int i = 0; i < particleNum; i++)
      {
        if (min_image(j, target.particles[i].r) < safeDistance && !(std::find(infected.begin(), infected.end(), i) != infected.end()))
        {
          infected.push_back(i);
          infectedColors.push_back(target.particles[i].strand);
        }
      }
    }
    inboxing();
    target.inboxing();
    std::vector<int> infectedStore = infected;
    cout << "Full infected list : "<<infected<<endl;
    cout << "N  "<< N << "   I   "<<infected.size()<<endl;
    if (infected.size() > N)
      std::cout << "Please reduce the safe distance multiplier. This feature is work in progress" << std::endl;

    for (int i=0;i<infected.size();i++){
      target.particles[infected[i]]+=particles[cluster[i]];
      cout << infected[i]<<",";
    }
    int p = infected.size();
    for (int j = 0; j < N; j++){
      for (int i = 0; i < particleNum; i++){
        if (min_image(j, target.particles[i].r) < safeDistance+0.2 && min_image(j, target.particles[i].r) > safeDistance && !(std::find(infected.begin(), infected.end(), i) != infected.end())){
          infected.push_back(i);
          infectedColors.push_back(target.particles[i].strand);
        }
      }
    }
    for (int i=p;i<N;i++){
      target.particles[infected[i]]+=particles[cluster[i]];
      cout << infected[i]<<",";
    }
    target.inboxing();
    target.customSeedForces(infectedStore);
    target.writeConfig();
    target.writeCrystalTopology();
    cout<<target.particles[418].r<<endl;
    return true;
  };

  // Computes the minimum distances for two particles with indices p and q
  template <typename A, typename B>double min_image(A p, B q){
    // valid only for cubic box.
    LR_vector p1, q1;
    if constexpr (is_same<A, LR_vector>::value)
    {
      p1 = p - npFloor(p / box.x) * box.x;
    }
    else if constexpr (is_same<A, int>::value)
    {
      p1 = particles[p].r - npFloor(particles[p].r / box.x) * box.x;
    }

    if constexpr (is_same<B, LR_vector>::value)
    {
      q1 = q - npFloor(q / box.x) * box.x;
    }
    else if constexpr (is_same<B, int>::value)
    {
      q1 = particles[q].r - npFloor(particles[q].r / box.x) * box.x;
    }

    p1 -= q1;
    p1 -= npRound(p1 / box.x) * box.x;
    return p1.module();
  }

  bool planeFitting(int *cluster, LR_vector center, LR_vector normal){

    return true;
  }

  bool writeCrystalTopology(std::string topology = ""){
    if (topology == "")
      topology = output + ".top";
    std::ofstream outputTop(topology);
    if (!outputTop.is_open())
      return false;
    outputTop << particleNum << " " << strands << std::endl;
    for (int i = 0; i < particleNum; i++)
    {
      outputTop << particles[i].strand << " ";
    }
    outputTop.close();
    return true;
  }

  bool writeConfig(std::string config = ""){

    if (config == "")
      config = output + ".dat";
    std::ofstream outputConfig(config);
    if (!outputConfig.is_open())
      return 21;
    outputConfig.precision(15);
    outputConfig << "t = 0" << std::endl;
    outputConfig << "b = " << box.x << " " << box.y << " " << box.z << std::endl;
    outputConfig << "E = 0 0 0" << std::endl;
    for (int i = 0; i < particleNum; i++)
    {
      outputConfig << particles[i].r.x << " ";
      outputConfig << particles[i].r.y << " ";
      outputConfig << particles[i].r.z << " ";
      outputConfig << particles[i].a1.x << " ";
      outputConfig << particles[i].a1.y << " ";
      outputConfig << particles[i].a1.z << " ";
      outputConfig << particles[i].a3.x << " ";
      outputConfig << particles[i].a3.y << " ";
      outputConfig << particles[i].a3.z << " ";
      outputConfig << "0 0 0 0 0 0" << std::endl;
    }
    return 0;
  }

private:
  string line, temp;
  istringstream ss;

  double subBoxing(double coordinate, double divisor){
    coordinate = std::remainder(coordinate, divisor);
    if (coordinate < 0)
      coordinate += divisor;
    return coordinate;
  };

  int readCrystalTopology(string topology){
    ifstream inputTop(topology);
    if (!inputTop.is_open())
      return 10;
    getline(inputTop, line);
    ss.clear();
    ss.str(line);
    ss >> particleNum;
    ss >> strands;
    particles.resize(particleNum);

    ss.clear();
    getline(inputTop, line);
    ss.str(line);

    for (i = 0; i < particleNum; i++)
    {
      ss >> temp;
      particles[i].id = i;
      particles[i].strand = stoi(temp);
    }
    return 0;
  }

  bool readDNAtopology(string topology){
    ifstream inputTop(topology);
    if (!inputTop.is_open())
      return false;
    getline(inputTop, line);
  }

  int readConfig(string config){
    ifstream inputConfig(config);
    if (!inputConfig.is_open())
      return 11;
    getline(inputConfig, line);
    getline(inputConfig, line);
    ss.clear();
    ss.str(line);
    ss >> temp;
    ss >> temp;
    ss >> box.x;
    ss >> box.y;
    ss >> box.z;
    getline(inputConfig, line);
    ss.clear();
    ss.str(line);
    ss >> temp;
    ss >> temp;
    ss >> energy.x;
    ss >> energy.y;
    ss >> energy.z;
    for (i = 0; i < particleNum; i++)
    {
      getline(inputConfig, line);
      ss.clear();
      ss.str(line);
      ss >> particles[i].r.x;
      ss >> particles[i].r.y;
      ss >> particles[i].r.z;
      ss >> particles[i].a1.x;
      ss >> particles[i].a1.y;
      ss >> particles[i].a1.z;
      ss >> particles[i].a3.x;
      ss >> particles[i].a3.y;
      ss >> particles[i].a3.z;
    }
    return 0;
  }

  int readPatches(string patches){
    ifstream inputPatches(patches);
    return 0;
  }

  int readParticles(string crystalpar){
    ifstream inputCrystal(crystalpar);
    return 0;
  }
};