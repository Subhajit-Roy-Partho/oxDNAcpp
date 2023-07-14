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

	// LR_vector LRremainder(LR_vector a, LR_vector b){
	// 	return (LR_vector){std::remainder(a.x,b.x),std::remainder(a.y,b.y),std::remainder(a.z,b.z)};
	// };

class Analysis {
public:
  int particleNum, strands, i;
  double safeMultiplier=2; // Multiplier with safe distance 
  string type,output;
  LR_vector box, energy;
  vector<Particle> particles;
  vector<vector<Traj>> traj; // For storing trajectory;
//   gsl_matrix *V = gsl_matrix_alloc(3, 3);
//   gsl_vector *S = gsl_vector_alloc(3);
//   gsl_vector *work = gsl_vector_alloc(3);
  Traj trajtemp;

  Analysis(string topology, string config, string type = "",
           string output = "output", string externalForces = "",
           string parameter1 = "", string parameter2 = "") {
    if (type == "crystal") {
      this->type = type;
      this->output=output;
      readCrystalTopology(topology);
      readConfig(config);
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
      particles[i].r.z = subBoxing(particles[i].r.x,box.z);
    }
  }

  bool pickAndPlace(int *cluster,int N,Analysis* target,LR_vector centralShift={0,0,0}) {
    Eigen::MatrixXd points(N,3);

    LR_vector center = CenterForIndex(cluster,N);
    double t = std::remainder((double)-7.1,(double)3);
    cout << "Redmainder = "<<t<<endl;
    // cout<<center<<endl;
    // cout <<particles[cluster[8]].r<<endl;
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
      // cout<<dist<<"\t";
      if(dist<safeDistance) safeDistance=dist;
    }
    // cout <<safeDistance<<endl;
    // cout<<npNorm((LR_vector) {points(0,0),points(0,1),points(0,2)} - (LR_vector){points(1,0),points(1,1),points(1,2)})<<endl;


    // for (int i=0; i<N;i++){
    //     particles[cluster[i]].r={points(i,0),points(i,1),points(i,2)};
    // }
    // writeCrystalTopology();
    // writeConfig();


    return true;
  };

  bool planeFitting(int *cluster, LR_vector center, LR_vector normal) {

    return true;
  }

    int writeCrystalTopology(string topology=""){
    if (topology=="") topology=output+".top";
    ofstream outputTop(topology);
    if(!outputTop.is_open()) return 20;
    outputTop<<particleNum<<" "<<strands<<endl;
    for(int i=0;i<particleNum;i++){
        outputTop<<particles[i].strand<<" ";
    }
    outputTop.close();
    return 0;
  }

  int writeConfig(string config=""){
    if (config=="") config=output+".dat";
    ofstream outputConfig(config);
    if(!outputConfig.is_open()) return 21;
    outputConfig.precision(15);
    outputConfig<<"t = 0"<<endl;
    outputConfig<<"b = "<<box.x<<" "<<box.y<<" "<<box.z<<endl;
    outputConfig<<"E = 0 0 0"<<endl;
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
        outputConfig<<"0 0 0 0 0 0"<<endl;
    }
    return 0;
  }

private:
  string line, temp;
  istringstream ss;

  double subBoxing(double coordinate,double divisor){
    coordinate = std::remainder(coordinate,divisor);
    if(coordinate<0) coordinate+=divisor;
    return 0; 
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