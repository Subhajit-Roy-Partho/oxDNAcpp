#include "Analysis.h"
#include <eigen3/Eigen/Dense>
#include <numeric>
#include <map>
using namespace std;


  bool Forces::add(string name, int particle, float stiff, LR_vector normal, LR_vector pos)
  {
    type.push_back(name);
    particles.push_back(particle);
    dir.push_back(normal);
    pos0.push_back(pos);
    this->stiff.push_back(stiff);
    return true;
  }
  bool Forces::addRepulsion(int particles, float stiff, LR_vector dir, float position)
  {
    type.push_back("repulsion_plane");
    this->particles.push_back(particles);
    this->stiff.push_back(stiff);
    this->dir.push_back(dir);
    this->position.push_back(position);
    return true;
  }
  bool Forces::addHarmonic(int particles, float stiff, float rate, LR_vector pos0, LR_vector dir)
  {
    type.push_back("trap");
    this->particles.push_back(particles);
    this->stiff.push_back(stiff);
    this->rate.push_back(rate);
    this->pos0.push_back(pos0);
    this->dir.push_back(dir);
    return true;
  }

  bool Forces::addHarmonicFromParticle(Particle *particle, float stiff, float rate)
  {
    return addHarmonic(particle->id, stiff, rate, particle->r);
  }

  bool Forces::save(string filename)
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
  };

  Analysis::Analysis(std::string topology, std::string config, std::string type, std::string output, std::string externalForces, std::string parameter1, std::string parameter2){
    this->type = type;
    this->output = output;
    this->topology = topology;
    this->config=config;
    
    if (type == "crystal"){
      readCrystalTopology(topology);
      readConfig(config);
      inboxing();
      // pickAndPlace();
    }
    if(type == "DNA"){
      if(!readDNAtopology(topology)){cout <<"Error reading DNA topology "<<endl;throw std::exception();}; 
      if(!readConfig(config)){cout <<"Error reading configuration file."<<endl;throw std::exception();};
      if(!shiftbox({0,0,0})){cout<<"Box overloaded from shift."<<endl;throw std::exception();};
    }

    if(type=="CCG"){
      if(!readCCGtopology(topology)){cout<<"Error reading CCG topolofy"<<endl;throw std::exception();};
      if(!readConfig(config)){cout<<"Error reading configuration file."<<endl;throw std::exception();};
    }
    // if(type == "newPSP"){
    // }
  }

  Analysis::~Analysis()
  {
    // gsl_matrix_free(V);
    // gsl_vector_free(S);
    // gsl_vector_free(work);
  }
  // Output the center for a number of index
  LR_vector Analysis::CenterForIndex(int *indexes, int N){
    LR_vector mean = {0, 0, 0};
    if(N==-1) return mean; // just ignore for default case.
    if(N==0){ // if N=0 return the center for the whole structure;
      for(int i=0;i<particleNum;i++){
        mean += particles[i].r;
      }
      mean /=particleNum; // This operation is only permitted if everything is in 1s quardant
      return mean;
    }
    for (int i = 0; i < N; i++){
      mean += particles[indexes[i]].r;
    }
    mean /= N;
    return mean;
  }

  LR_vector Analysis::CenterForIndex(vector<int> indexes){
    LR_vector mean = {0, 0, 0};
    int N = indexes.size();
    for (int i = 0; i < N; i++){
      mean += particles[indexes[i]].r;
    }
    mean /= N;
    return mean;
  }

LR_vector Analysis::CenterForIndex(int N){
  LR_vector mean = {0,0,0};
  for(i=0;i<particleNum;i++){
    mean += particles[i].r;
  }
  mean /= particleNum;
  if(N==-1){//mean of the whole structure
    return mean;
  }
  mean -= particles[i].r; // module of this would be the distance of the particle from the center.
  return mean;
}

  LR_vector Analysis::NormalToPlane(Eigen::MatrixXd points)
  {
    Eigen::Matrix3d m(3, 3);
    m << points.transpose() * points;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
    m = svd.matrixV();
    return (LR_vector){m(0, 2), m(1, 2), m(2, 1)};
  };
  void Analysis::inboxing(LR_vector center)
  {
    for (int i = 0; i < particleNum; i++)
    {
      particles[i].r.x = subBoxing(particles[i].r.x, box.x);
      particles[i].r.y = subBoxing(particles[i].r.y, box.y);
      particles[i].r.z = subBoxing(particles[i].r.z, box.z);
    }
  }

  bool Analysis::customSeedForces(std::vector<int> ids)
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

  bool Analysis::correctA(string newFile){
    Analysis newPar(topology,newFile,"crystal");
    for(int i=0;i<particleNum;i++){
      particles[i].a1=newPar.particles[i].a1;
      particles[i].a3=newPar.particles[i].a3;
    }
    return true;
  }

  bool Analysis::randomReplaceColor(int originalColor, int newColor, vector<int> ignore, int N)
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

  bool Analysis::randomReplacePosition(int originalColor, Particle *newParticle, vector<int> ignore)
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
  bool Analysis::pickAndPlace(int *cluster, int N, string tname, LR_vector centralShift)
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
    // cout << "Full infected list : "<<infected<<endl;
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
  template <typename A, typename B>double Analysis::min_image(A p, B q){
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

  bool Analysis::planeFitting(int *cluster, LR_vector center, LR_vector normal){

    return true;
  }

  bool Analysis::writeCrystalTopology(std::string topology){
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

  bool Analysis::writeDNAtopology(std::string topology){
    return true;
  }

  bool Analysis::writeConfig(std::string config){

    if (config == "")
      config = output + ".dat";
    std::ofstream outputConfig(config,ios::trunc);
    if (!outputConfig.is_open())
      return false;
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
    outputConfig.close();
    return true;
  }

  bool Analysis::writeMGL(std::string config){
    if(config == "") config = output+".mgl";
    std::ofstream outputMGL(config,ios::trunc);
    if(!outputMGL.is_open()) return false;
    outputMGL<<".Box:"<<box<<std::endl;
    for(i=0;i<strands;i++){
      auto part = particles[i*particlePerStrand+particlePerStrand-1].r;
      outputMGL<<part.x<<" "<<part.y<<" "<<part.z<<" @ 18 C[]"<<std::endl;
      for(int j=particlePerStrand-2;j>=0;j--){
      }
    }
    return true;
  }

  double Analysis::subBoxing(double coordinate, double divisor){
    coordinate = std::remainder(coordinate, divisor);
    if (coordinate < 0) coordinate += divisor;
    return coordinate;
  };

  bool Analysis::readCrystalTopology(string topology){
    ifstream inputTop(topology);
    if (!inputTop.is_open())
      return false;
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
    return true;
  }

  bool Analysis::readDNAtopology(string topology){
    ifstream inputTop(topology);
    if (!inputTop.is_open())
      return false;
    getline(inputTop, line);
    ss.clear();
    ss.str(line);
    ss >> particleNum;
    ss >> strands;
    particles.resize(particleNum);
    for(i=0;i<particleNum;i++){
      ss.clear();
      getline(inputTop,line);
      ss.str(line);
      ss>>particles[i].strand;
      ss>>particles[i].name;
      ss>>tempInt;particles[i].connector.push_back(tempInt);
      ss>>tempInt;particles[i].connector.push_back(tempInt);
    }
    return true;
  }

  bool Analysis::readCCGtopology(string topology){
    ifstream inputTop(topology);
    if (!inputTop.is_open())
      return false;
    getline(inputTop, line);
    ss.clear();
    ss.str(line);
    ss >> particleNum;
    ss >> strands;
    particles.resize(particleNum);
    particlePerStrand=particleNum/strands;
    for(i=0;i<particleNum;i++){
      ss.clear();
      getline(inputTop,line);
      ss.str(line);
      ss>>temp;
      ss>>particles[i].strand;
      ss>>particles[i].color;
      ss>>particles[i].radius;
    }
    return true;
  }

  bool Analysis::readConfig(string config){
    ifstream inputConfig(config);
    if (!inputConfig.is_open())
      return false;
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

    int count=0;
    if(particleNum==0){
      particles.resize(99999);
      while(getline(inputConfig, line)){
        if(line.empty()||line[0]=='#') continue;
        ss.clear();
        ss.str(line);
        ss >> particles[count].r.x;
        ss >> particles[count].r.y;
        ss >> particles[count].r.z;
        ss >> particles[count].a1.x;
        ss >> particles[count].a1.y;
        ss >> particles[count].a1.z;
        ss >> particles[count].a3.x;
        ss >> particles[count].a3.y;
        ss >> particles[count].a3.z;
        count++;
      }
      particles.resize(count);
      particleNum=count;
      return true;
    }

    for (i = 0; i < particleNum; i++)
    {
      getline(inputConfig, line);
      // if(line.empty()||line[0]=='#') continue;
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
    return true;
  }

  bool Analysis::readPatches(string patches){
    ifstream inputPatches(patches);
    return false;
  }

  bool Analysis::readParticles(string crystalpar){
    ifstream inputCrystal(crystalpar);
    return false;
  }

bool Analysis::shiftbox(LR_vector shift){
  if(shift==LR_vector({0,0,0})){
    LR_vector minimum={0,0,0};
    for(i=0;i<particleNum;i++){
      if(minimum.x>particles[i].r.x) minimum.x= particles[i].r.x;
      if(minimum.y>particles[i].r.y) minimum.y= particles[i].r.y;
      if(minimum.z>particles[i].r.z) minimum.z= particles[i].r.z;
    }
    shiftbox(minimum.abs());
    
  }else{
    for(i=0;i<particleNum;i++){
      particles[i].r+=shift;
    }
    return false;
  }
  return true;
}

bool Analysis::testBoxOverloaded(){
  LR_vector minimum ={0,0,0},maximum = box;
  for(i=0;i<particleNum;i++){
    if(minimum.x>particles[i].r.x) minimum.x= particles[i].r.x;
    if(maximum.x<particles[i].r.x) return false;
    if(minimum.y>particles[i].r.y) minimum.y= particles[i].r.y;
    if(maximum.y<particles[i].r.y) return false;
    if(minimum.z>particles[i].r.z) minimum.z= particles[i].r.z;
    if(maximum.z<particles[i].r.x) return false;
  }
  if (minimum.x<0) return false;
  if (minimum.y<0) return false;
  if (minimum.z<0) return false;
  return true;
}

bool Analysis::generatePSP(Analysis *PSP,vector<vector<int>> ids,vector<int> colors,vector<double> radius,int numNeighbour,double fixedSpringConstant){
  if(colors.size()<ids.size() || colors.size()>ids.size()+1){
    cout << "Invalid number of colors are passed."<<endl;
    return false;
  }
  // Setting up all the parameters of the particle
  PSP->strands=1;
  PSP->particleNum=ids.size()+1;
  PSP->particles.resize(PSP->particleNum);
  PSP->particlePerStrand=ids.size()+1;
  PSP->box={50,50,50};

  vector<double> v; // stores the distance between the clusters
  int j=0;
  PSP->particles[ids.size()].r = CenterForIndex((int)-1);
  for(int j=0;j<ids.size();j++){
    v.clear();
    PSP->particles[j].color=colors[j];
    PSP->particles[j].radius=radius.size()==2?radius[0]:radius[j];


    PSP->particles[j].r=CenterForIndex(ids[j]);
    for(i=0;i<ids.size();i++)  v.push_back(((PSP->particles[j].r-CenterForIndex(ids[i]))).module());
    auto sorted = sort_indexes(v);// sort the index 
    
    for(int p=0;p<numNeighbour+1;p++){ // connect only numNeighbour particles
      if(sorted[p]==j) continue;
      PSP->particles[j].connector.push_back(sorted[p]);
      PSP->particles[j].eqRadius.push_back(v[sorted[p]]);
      if(fixedSpringConstant!=0) PSP->particles[j].spring.push_back(fixedSpringConstant);// spring is hardcoded for 100
    }
    PSP->particles[j].connector.push_back(ids.size()); // Add the last central particle to the list
    PSP->particles[j].eqRadius.push_back((PSP->particles[ids.size()].r-PSP->particles[j].r).module());
    if(fixedSpringConstant!=0) PSP->particles[j].spring.push_back(fixedSpringConstant); // spring is hardcoded for 100
    // cout<<PSP->particles[j].eqRadius<<endl;
  }
  PSP->particles[ids.size()].color=colors.size()== ids.size()?100:colors[ids.size()];
  PSP->particles[ids.size()].radius=radius.size()==2?radius[1]:radius[ids.size()];

  int max=0;
  for(int j=0;j<3;j++){   
    for(i=0;i<PSP->particleNum;i++){
      if(std::ceil(PSP->particles[i].r[j])>max) max=std::ceil(PSP->particles[i].r[j]);
    }
    PSP->box[j]=max+2;
  }
  return true;
}


template <typename T> vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

bool Analysis::writeCCGtopology(string topology){
  if (topology == "") topology = output + ".top";
  ofstream outputTop(topology,ios::trunc);
  if(!outputTop.is_open()) return false;

  outputTop.precision(15);
  outputTop<<particleNum<<" "<<strands<<" "<<particleNum<<" 0 0 0\n";
  for (i=0;i<particleNum;i++){
    if(particles[i].connector.size()==0){
      outputTop<<"-2 "<<particles[i].strand<<" "<<particles[i].color<<" "<<particles[i].radius;
    }else{
      outputTop<<"-2 "<<particles[i].strand<<" "<<particles[i].color<<" "<<particles[i].radius<<" ";
    }
    for(int j=0;j<particles[i].connector.size();j++){
      if(j==particles[i].connector.size()-1){outputTop<<particles[i].connector[j]<<" "<<particles[i].spring[j]<<" "<<particles[i].eqRadius[j];}
      else{outputTop<<particles[i].connector[j]<<" "<<particles[i].spring[j]<<" "<<particles[i].eqRadius[j]<<" ";};
    }
    outputTop<<std::endl;
  }
  outputTop.close();
  return true;
}

bool Analysis::writeCCGviewTopology(string topology){
  if(topology=="") topology = "view"+output+".top";
  ofstream outputTop(topology,ios::trunc);
  if(!outputTop.is_open()) return false;
  outputTop.precision(0);
  outputTop<<particleNum<<" "<<strands<<endl;

  std::map<int, std::string> fakeNucluo = {{0,"A"},
                                            {1,"T"},
                                            {2,"G"},
                                            {3,"C"}};
  bool first=true;
  for(i=0;i<particleNum;i++){
    if(first){
      outputTop<<particles[i].strand+1<<" "<<fakeNucluo[particles[i].strand%4]<<" "<<-1<<" "<<i+1<<endl;
      first=false;
    }else if(particles[i].connector.size()==0){
      outputTop<<particles[i].strand+1<<" "<<fakeNucluo[particles[i].strand%4]<<" "<<i-1<<" "<<-1<<endl;
      first =true;
    }else outputTop<<particles[i].strand+1<<" "<<fakeNucluo[particles[i].strand%4]<<" "<<i-1<<" "<<i+1<<endl;
  }
  outputTop.close();
  return true;
}

bool Analysis::populate(int num,double seperator){
  int totPar=particleNum;
  LR_vector minSize = box-2+seperator;
  // Genral parameters
  strands*=num;
  particleNum*=num;
  particles.resize(particleNum);
  
  int dim=std::ceil(std::pow(num,1.0/3.0));
  box=minSize*dim+2;

  LR_vector translate={0,0,0};
  int currentNum=1;
  buildingParticles:
    for(i=0;i<dim;i++){
      for(int j=0;j<dim;j++){
        for(int k=0;k<dim;k++){
          if(i==0 && j==0 && k==0) continue;
          LR_vector shift={minSize.x*k,minSize.y*j,minSize.z*i};
          for(int p=0;p<totPar;p++){
            int index=p+k*totPar+j*dim*totPar+i*dim*dim*totPar;
            particles[index] =particles[p];
            particles[index].r+=shift;
            particles[index].strand=currentNum;
            // std::for_each(particles[index].connector.begin(),particles[index].connector.end(),[](& d){d+=fcurrentNum*13;});
            for(int m=0;m<particles[index].connector.size();m++) particles[index].connector[m]+=currentNum*particlePerStrand;
          }
          currentNum+=1;
          if(currentNum==num) goto exit;
        }
      }
    }
    return true;
  exit:
    return true;

  return true;
}

bool Analysis::boxToCubic(){
  double max=0;
  for(i=0;i<3;i++){
    if(max<box[i])max=box[i];
  }
  box=(LR_vector){max,max,max};
  return true;
}
