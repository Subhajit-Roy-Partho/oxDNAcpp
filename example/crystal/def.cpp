/* ############## Errors #####################
10 - error opening input Topology
11 - error opening input Config
*/

#include "LR_vector.h"
#include <cmath>
#include <fstream>
// #include <gsl/gsl_linalg.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <eigen3/Eigen/Dense>


using namespace std;

template <typename S>
std::ostream& std::operator<<(std::ostream& os,
                    const std::vector<S>& vector)
{
    // Printing all the elements
    // using <<
    for (auto element : vector) {
        os << element << " ";
    }
    return os;
}

class Particle {
public:
  int id, color, strand;
  LR_vector r = {0, 0, 0}, a1, a3;
};

struct Traj {
  LR_vector r;
  LR_vector a1;
  LR_vector a3;
};

// Tried to create numpy like function valid for atleast LR_vector
template <typename A, std::size_t N> 
A npMean(A (&vector)[N]) {
  A result = vector[0];
  for (int i = 1; i < N; i++) {
    result += vector[i];
  }
  return result / N;
};

template <typename A>
double npNorm(A vector, int norm =2){
  if(is_same<A,LR_vector>::value)
    return pow(pow(vector.x,norm) + pow(vector.y,norm)+pow(vector.z,norm),1.0/(double)norm);
  
  return 0;
}

template <typename A>
A npFloor(A vector){
  if(is_same<A,LR_vector>::value) return (LR_vector){std::floor(vector.x),std::floor(vector.y),std::floor(vector.z)};
}

template <typename A>
A npRound(A vector){
  if(is_same<A,LR_vector>::value) return (LR_vector){std::round(vector.x),std::round(vector.y),std::round(vector.z)};
}

class Analysis {
public:
  int particleNum, strands, i;
  double safeMultiplier=1.4; // Multiplier with safe distance 
  std::string type,output;
  LR_vector box, energy;
  std::vector<Particle> particles;
  std::vector<std::vector<Traj>> traj; // For storing trajectory;
//   gsl_matrix *V = gsl_matrix_alloc(3, 3);
//   gsl_vector *S = gsl_vector_alloc(3);
//   gsl_vector *work = gsl_vector_alloc(3);
  Traj trajtemp;

  Analysis(std::string topology, std::string config, std::string type = "",
           std::string output = "output", std::string externalForces = "",
           std::string parameter1 = "", std::string parameter2 = "") {
    if (type == "crystal") {
      this->type = type;
      this->output=output;
      readCrystalTopology(topology);
      readConfig(config);
      inboxing();
      // pickAndPlace();
    }
  }
  ~Analysis() {
    // gsl_matrix_free(V);
    // gsl_vector_free(S);
    // gsl_vector_free(work);
  }
  // Output the center for a number of index
  LR_vector CenterForIndex(int *indexes,int N) {
    LR_vector mean = {0, 0, 0};
    for (int i = 0; i < N; i++) {
      mean += particles[indexes[i]].r;
    }
    mean /= N;
    return mean;
  }

  LR_vector NormalToPlane(Eigen::MatrixXd points){
    Eigen::Matrix3d m(3,3);
    m << points.transpose()*points;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(m,Eigen::ComputeFullU| Eigen::ComputeFullV);
    m=svd.matrixV();
    return (LR_vector) {m(0,2),m(1,2),m(2,1)};

  };
  void inboxing(LR_vector center ={0,0,0}){
    for(int i=0;i<particleNum;i++){
      particles[i].r.x = subBoxing(particles[i].r.x,box.x);
      particles[i].r.y = subBoxing(particles[i].r.y,box.y);
      particles[i].r.z = subBoxing(particles[i].r.z,box.z);
    }
  }

  bool pickAndPlace(int *cluster,int N,Analysis* target,LR_vector centralShift={0,0,0}) {
    Eigen::MatrixXd points(N,3);

    LR_vector center = CenterForIndex(cluster,N);
    double t = std::remainder((double)-7.1,(double)3);
    for(int i=0;i<N;i++){
      LR_vector some = particles[cluster[i]].r -particles[cluster[i+1]].r;
      // cout<<some.module()<<endl;
      particles[cluster[i]].r-=center+centralShift;
      points(i,0)=particles[cluster[i]].r.x;
      points(i,1)=particles[cluster[i]].r.y;
      points(i,2)=particles[cluster[i]].r.z;
    }

    LR_vector normal = NormalToPlane(points);
    // cout<< normal<<endl;

    Eigen::Matrix3d Rot;
    Rot=Eigen::AngleAxisd(acos(normal.x)-M_PI,Eigen::Vector3d::UnitY())*
        Eigen::AngleAxisd(M_PI-acos(normal.y),Eigen::Vector3d::UnitX());

    // cout<<(LR_vector) {points(0,0),points(0,1),points(0,2)}<<endl;
    points= points*Rot;
    double safeDistance=99999999999999999;
    double dist;
    for (int i=0;i<N-1;i++){
      dist = npNorm((LR_vector) {points(i,0),points(i,1),points(i,2)} - (LR_vector){points(i+1,0),points(i+1,1),points(i+1,2)});
      if(dist<safeDistance) safeDistance=dist;
    }
    safeDistance *=safeMultiplier;
    // std::cout <<safeDistance<<std::endl;
    std::vector <int> infected,infectedColors;
    // cout <<min_image(1,particles[2].r)<<endl;
    for(int j=0;j<N;j++){
      for(int i=0;i<particleNum;i++){
        if(min_image(j,target->particles[i].r)<safeDistance && !(std::find(infected.begin(),infected.end(),i)!=infected.end())) {
          infected.push_back(i);
          infectedColors.push_back(target->particles[i].strand);
        }
      }
    }
    if(infected.size()>N) std::cout<<"Please reduce the safe distance multiplier. This feature is work in progress"<<std::endl;

    for(int i=0;i<N;i++){
      int tempColor=particles[cluster[i]].strand;
      std::cout << "Color = "<<tempColor<<std::endl;
      std::vector<int>::iterator itr = std::find(infectedColors.begin(),infectedColors.end(),tempColor);
      if(itr!=infectedColors.cend()){
        int index = std::distance(infectedColors.begin(),itr);
        cout<< "Index = "<<index<<endl;
        cout<<"Prev:\n"<<infectedColors<<endl;
        infectedColors.erase(infectedColors.begin()+index);
        infected.erase(infected.begin()+index);
        cout <<"Now:\n"<<infectedColors<<endl;
      }else{
        cout<<"Outside"<<endl;
      }
    }
    cout <<infected<<endl;

    // writeConfig();
    return true;
  };

  // Computes the minimum distances for two particles with indices p and q
  template <typename A,typename B>
  double min_image(A p, B q){
    //valid only for cubic box.
    LR_vector p1,q1;
    if constexpr(is_same<A,LR_vector>::value){
      p1=p-npFloor(p/box.x)*box.x;
    }else if constexpr(is_same<A,int>::value){
      p1 = particles[p].r-npFloor(particles[p].r/box.x)*box.x;
    }
    
    if constexpr(is_same<B,LR_vector>::value){
      q1=q-npFloor(q/box.x)*box.x;
    }else if constexpr(is_same<B,int>::value){
      q1 = particles[q].r-npFloor(particles[q].r/box.x)*box.x;
    }

    p1-=q1;
    p1 -= npRound(p1/box.x)*box.x;
    return p1.module();
  }

  bool planeFitting(int *cluster, LR_vector center, LR_vector normal) {

    return true;
  }

    int writeCrystalTopology(std::string topology=""){
    if (topology=="") topology=output+".top";
    std::ofstream outputTop(topology);
    if(!outputTop.is_open()) return 20;
    outputTop<<particleNum<<" "<<strands<<std::endl;
    for(int i=0;i<particleNum;i++){
        outputTop<<particles[i].strand<<" ";
    }
    outputTop.close();
    return 0;
  }

  int writeConfig(std::string config=""){
    if (config=="") config=output+".dat";
    std::ofstream outputConfig(config);
    if(!outputConfig.is_open()) return 21;
    outputConfig.precision(15);
    outputConfig<<"t = 0"<<std::endl;
    outputConfig<<"b = "<<box.x<<" "<<box.y<<" "<<box.z<<std::endl;
    outputConfig<<"E = 0 0 0"<<std::endl;
    for (int i=0;i<particleNum;i++){
        outputConfig<<particles[i].r.x<<" ";
        outputConfig<<particles[i].r.y<<" ";
        outputConfig<<particles[i].r.z<<" ";
        outputConfig<<particles[i].a1.x<<" ";
        outputConfig<<particles[i].a1.y<<" ";
        outputConfig<<particles[i].a1.z<<" ";
        outputConfig<<particles[i].a3.x<<" ";
        outputConfig<<particles[i].a3.y<<" ";
        outputConfig<<particles[i].a3.z<<" ";
        outputConfig<<"0 0 0 0 0 0"<<std::endl;
    }
    return 0;
  }

private:
  string line, temp;
  istringstream ss;

  double subBoxing(double coordinate,double divisor){
    coordinate = std::remainder(coordinate,divisor);
    if(coordinate<0) coordinate+=divisor;
    return coordinate; 
  };

  int readCrystalTopology(string topology) {
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

    for (i = 0; i < particleNum; i++) {
      ss >> temp;
      particles[i].id = i;
      particles[i].strand = stoi(temp);
    }
    return 0;
  }

  int readConfig(string config) {
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
    for (i = 0; i < particleNum; i++) {
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

  int readPatches(string patches) {
    ifstream inputPatches(patches);
    return 0;
  }

  int readParticles(string crystalpar) {
    ifstream inputCrystal(crystalpar);
    return 0;
  }
};



// template <typename A>

// Selected Points
// 1070,134,145,155,1004,1385,560,1291,904,136,1112,126,1355,1850,1384,686,392,1605,618,1024,2027,474,1640,1720,1831,387,1504,947,1832,1262,676,1067,195,121,734,245,217


    // double points[] = {158720.15575206, 42724.03921793,  56622.47200362,
    //                    42724.03921793,  132381.4182789,  -83288.45034046,
    //                    56622.47200362,  -83288.45034046, 100855.62941231};
    // gsl_matrix_view U = gsl_matrix_view_array(points, 3, 3);
    // int pass = gsl_linalg_SV_decomp(&U.matrix, V, S, work);
    // gsl_matrix_fprintf(stdout, V, "%g");