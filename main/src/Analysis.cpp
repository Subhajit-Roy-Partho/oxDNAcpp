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
      if(!readCrystalTopology(topology)){cout <<"Error reading Patchy Crystal topology "<<endl;throw std::exception();};
      if(!readConfig(config)){cout <<"Error reading configuration file."<<endl;throw std::exception();};
      // if(!readCrystalPatches())
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

    if(type=="trajPHB"){
      if(!readPHBTopology(topology)){cout<<"Error reading PHB topology"<<endl;throw std::exception();};
      if(!readTrajectory(config)){cout<<"Error reading trajectory file."<<endl;throw std::exception();};
    }

    size_t pos=0;
    if((pos=type.find("new"))!=std::string::npos){ //make sure if new keyword is attached keep the type same and don't do extra stuff
      type.erase(0,pos+3);
      this->type=type;
    }
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

  bool Analysis::writeMGL(std::string config,double centralRadius,double patchRadius,bool modern){
    if(config == "") config = output+".mgl";
    std::ofstream outputMGL(config,ios::trunc);
    if(!outputMGL.is_open()) return false;
    outputMGL.precision(15);
    outputMGL<<".Box:"<<box<<std::endl;
    for(i=0;i<particleNum;i++){
      outputMGL<<particles[i].r.x<<" "<<particles[i].r.y<<" "<<particles[i].r.z<<" @ ";
      if(centralRadius==0){
        outputMGL<<particles[i].radius<<" ";
      }else{
        outputMGL<<centralRadius<<" ";
      };
      outputMGL<<"C["<<centralColorMap.find(particles[i].strand)->second<<"] M ";
      for(int j=0;j<particles[i].patches.size();j++){
        auto temp=sourcePatch[particles[i].patches[j]];
        if(patchRadius==0) patchRadius=temp.radius;
        outputMGL<<temp.position.x<<" "<<temp.position.y<<" "<<temp.position.z<<" "<<patchRadius<<" C["<<colorMap.find(temp.color)->second<<"] ";
      }
      outputMGL<<std::endl;
    }
    return true;
  }

  bool Analysis::writeMGLtraj(std::string topology, int start, int end, int step,bool truncate,double patchRadius){
    // Get end if end is -1
    if(end==-1) end=traj.size();
    if(topology=="") topology = output+".mgl";

    std::ofstream outputMGL(topology,ios::app);
    if(!outputMGL.is_open()) return false;

    if(truncate){
      std::ofstream outputMGL(topology,ios::trunc);
      if(!outputMGL.is_open()) return false;
    }

    //get unique colors
    std::vector<int> uniqueColors;
    for(uint p=0;p<sourcePatch.size();p++){
      uniqueColors.push_back(sourcePatch[p].color);
    }
    std::set<int> unique(uniqueColors.begin(),uniqueColors.end());
    uniqueColors.assign(unique.begin(),unique.end());
    // for(int value:uniqueColors){
    //   cout<<value<<endl;
    // }

    auto colorRGB = generateRandomColors(uniqueColors.size());
    auto colorParticle = generateRandomColors(patchConfig.size());

    outputMGL.precision(15);
    for(i=start;i<end;i+=step){
      outputMGL<<".Box:"<<box<<std::endl;
      for(int j=0;j<particleNum;j++){
        outputMGL<<traj[i].r[j].x<<" "<<traj[i].r[j].y<<" "<<traj[i].r[j].z<<" @ ";
        outputMGL<<particles[j].radius<<" ";
        if(particles[j].color==100){
          outputMGL<<"C[255,255,255,1] ";
          outputMGL<<std::endl;
          continue;
        }
        LR_vector particleColor = colorParticle[particles[j].color];
        outputMGL<<"C["<<(int)particleColor.x<<","<<(int)particleColor.y<<","<<(int)particleColor.z<<",1] M ";
        for(int k=0;k<particles[j].patches.size();k++){
          auto temp=sourcePatch[particles[j].patches[k]];
          auto it = std::find(uniqueColors.begin(),uniqueColors.end(),temp.color);
          if(it!=uniqueColors.end()){
            int index = std::distance(uniqueColors.begin(),it);
            LR_vector patchColor = colorRGB[index];
            outputMGL<<temp.position.x<<" "<<temp.position.y<<" "<<temp.position.z<<" "<<patchRadius<<" C["<<(int)patchColor.x<<","<<(int)patchColor.y<<","<<(int)patchColor.z<<",1] ";
          }
        }
        outputMGL<<std::endl;
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

  bool Analysis::readCrystalPatches(std::string patch){
    ifstream inputPatch(patch);
    if(!inputPatch.is_open()) return false;
    bool open=false;
    int count = 0;
    // std::string delimiter = "=";
    Patch tempPatch;
    while(getline(inputPatch,line)){
      if(!line.size()) continue; //ignore empty lines 
      if(line[0]=='#') continue; //ignore comments 
      if(line.find("{") != std::string::npos){
        open=true;
        continue;
      };
      if(line.find("}") != std::string::npos){
        if(open){
          sourcePatch.push_back(tempPatch);
          count++;
        }
        open=false;
        continue;
      }
      auto splitted = npSplit(line,"=");
      if(splitted[0]=="id") if(stoi(splitted[1])!=count) return false;
      if(splitted[0]=="color") tempPatch.color=stoi(splitted[1]);
      if(splitted[0]=="strength") tempPatch.strength=stoi(splitted[1]);
      if(splitted[0]=="position"){
        splitted = npSplit(splitted[1],",");
        tempPatch.position=(LR_vector){stod(splitted[0]),stod(splitted[1]),stod(splitted[2])};
        // cout<<tempPatch.position<<endl;
      }
      if(splitted[0]=="a1"){
        splitted=npSplit(splitted[1],",");
        tempPatch.a1=(LR_vector){stod(splitted[0]),stod(splitted[1]),stod(splitted[2])};
      }
      if(splitted[0]=="a2"){
        splitted=npSplit(splitted[1],",");
        tempPatch.a2=(LR_vector){stod(splitted[0]),stod(splitted[1]),stod(splitted[2])};
      }
    }
    return true;
  };

  bool Analysis::readCrystalParticlePatchyConfig(std::string file){
    ifstream inputConfig(file);
    if(!inputConfig.is_open()) return false;
    bool open=false;
    int count = 0;
    std::vector<int> tempConfig;
    while(getline(inputConfig,line)){
      if(!line.size()) continue; //ignore empty lines 
      if(line[0]=='#') continue; //ignore comments 
      if(line.find("{") != std::string::npos){
        open=true;
        continue;
      };
      if(line.find("}") != std::string::npos){
        if(open){
          patchConfig.push_back(tempConfig);
          tempConfig.clear();
          count++;
        }
        open=false;
        continue;
      }
      auto splitted = npSplit(line,"=");
      if(splitted[0]=="type") if(stoi(splitted[1])!=count) return false;
      if(splitted[0]=="patches"){
        splitted = npSplit(splitted[1],",");
        for(int i=0;i<splitted.size();i++){
          tempConfig.push_back(stoi(splitted[i]));
        }
        for(i=0;i<particleNum;i++)
          if(particles[i].strand==count) particles[i].patches=tempConfig;
      }
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

  bool Analysis::readPHBTopology(string topology){
    ifstream inputTop(topology);
    if (!inputTop.is_open())
      return false;
    getline(inputTop, line);
    ss.clear();
    ss.str(line);
    ss >> particleNum;
    ss >> strands;
    particles.resize(particleNum);
    particlePerStrand=-1;
    Patch tempPatch;
    int i=0;
    while(std::getline(inputTop,line)){
      // if(i>=particleNum-1) break;
      std::stringstream body(line);
      if(line.empty()||line[0]=='#') continue;
      if(line[0]=='i'){
        if(line[1]=='P'){
          int j=0;
          while(body.tellg()!=-1){
            body>>temp;
            if(j==2) tempPatch.color=std::stoi(temp);
            if(j==3) tempPatch.strength=std::stod(temp);
            if(j==4) tempPatch.position.x=std::stod(temp);
            if(j==5) tempPatch.position.y=std::stod(temp);
            if(j==6) tempPatch.position.z=std::stod(temp);
            if(j==7) tempPatch.a1.x=std::stod(temp);
            if(j==8) tempPatch.a1.y=std::stod(temp);
            if(j==9) tempPatch.a1.z=std::stod(temp);
            if(j==10) tempPatch.a2.x=std::stod(temp);
            if(j==11) tempPatch.a2.y=std::stod(temp);
            if(j==12) tempPatch.a2.z=std::stod(temp);
            j++;
          }
          sourcePatch.push_back(tempPatch);
          continue;
        }
        if(line[1]=='C'){
          std::vector<int> tempConfig;
          int j=0;
          while(body.tellg()!=-1){
            body>>temp;
            if(j>1) tempConfig.push_back(std::stoi(temp));
            j++;
          }
          patchConfig.push_back(tempConfig);
          continue;
        }
      }
      // cout<<line<<endl;
      int j=0;
      // std::stringstream body(line);
      while(body.tellg()!=-1){
        body>>temp;
        if(j==0) particles[i].type=std::stoi(temp);
        if(j==1) particles[i].strand=std::stoi(temp);
        if(j==2) {
          particles[i].color=std::stoi(temp);
          if(particles[i].color!=100)
          particles[i].patches=patchConfig[particles[i].color];
        };
        if(j==3) particles[i].radius=std::stod(temp);
        j++;
      }
      i++;
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

  bool Analysis::readTrajectory(std::string trajectory,unsigned int start,int end, unsigned int step){
    ifstream inputTraj(trajectory);
    if(!inputTraj.is_open()) return false;

    traj.clear();
    temptraj.updateParticleNumber(particleNum);
    if(end==-1) while(getline(inputTraj,line)){
      ss.clear();ss.str(line);ss>>temp;ss>>temp;ss>>temp;
      temptraj.time=stoi(temp);
      getline(inputTraj,line); //skip the box line
      if(box==LR_vector({0,0,0})){
        ss.clear();ss.str(line);
        ss>>temp;ss>>temp;
        ss>>box.x;ss>>box.y;ss>>box.z;
      }
      getline(inputTraj,line); // skip the energy line
      for(i=0;i<particleNum;i++){
        getline(inputTraj,line);
        ss.clear();ss.str(line);
        ss>>temptraj.r[i].x;
        ss>>temptraj.r[i].y;
        ss>>temptraj.r[i].z;
        ss>>temptraj.a1[i].x;
        ss>>temptraj.a1[i].y;
        ss>>temptraj.a1[i].z;
        ss>>temptraj.a3[i].x;
        ss>>temptraj.a3[i].y;
        ss>>temptraj.a3[i].z;
      }
      traj.push_back(temptraj);
      if(step>1) skipLines(inputTraj,(step-1)*(particleNum+3));
    }else{
      traj.resize((end-start)/step);
      traj.shrink_to_fit();
      for(i=start;i<end;i+=step){
        getline(inputTraj,line); // time
        temptraj.time=stoi(line);
        getline(inputTraj,line); //skip the box line
        if(box==LR_vector({0,0,0})){
          ss.clear();ss.str(line);
          ss>>temp;ss>>temp;
          ss>>box.x;ss>>box.y;ss>>box.z;
        }
        getline(inputTraj,line); // skip the energy line
        for(int j=0;j<particleNum;j++){
          getline(inputTraj,line);
          ss.clear();ss.str(line);
          ss>>temptraj.r[j].x;
          ss>>temptraj.r[j].y;
          ss>>temptraj.r[j].z;
          ss>>temptraj.a1[j].x;
          ss>>temptraj.a1[j].y;
          ss>>temptraj.a1[j].z;
          ss>>temptraj.a3[j].x;
          ss>>temptraj.a3[j].y;
          ss>>temptraj.a3[j].z;
        }
        traj.push_back(temptraj);
      }
    }
    inputTraj.close();
    return true;
  }

  // bool Analysis::readPatches(string patches){
  //   ifstream inputPatches(patches);
  //   return false;
  // }

  // bool Analysis::readParticles(string crystalpar){
  //   ifstream inputCrystal(crystalpar);
  //   return false;
  // }

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

bool Analysis::generatePSP(Analysis *PSP,vector<vector<int>> ids,vector<int> &colors,vector<double> &radius,int numNeighbour,double fixedSpringConstant){
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

bool Analysis::addNewType(LR_vector shift,vector<int> colors,vector<double> radius){
  if(colors.size()==particlePerStrand-1){// As the central particle is always 100 so prevent redundency
    colors.push_back(100);
  } else if(colors.size()<particlePerStrand){// if less number of colors are provided the process will fails
    cout << "Invalid number of colors, operation failed"<<endl;
    return false;
  }
  strands+=1;
  particleTypes+=1;
  particleNum+=particlePerStrand;
  particles.resize(particleNum);
  for(i=0;i<3;i++) box[i]=box[i]*particleTypes/(particleTypes-1); //go over x,y,z of the box
  LR_vector minSize=(box/particleTypes);
  // cout<<minSize.multiplyEach(shift)<<endl;
  for(i=0;i<particlePerStrand;i++){
    int index = (strands-1)*particlePerStrand+i;
    particles[index]=particles[i];
    particles[index].color=colors[i];
    particles[index].r+=minSize.multiplyEach(shift*1.3);
    particles[index].strand=strands;
    for(int m=0;m<particles[index].connector.size();m++) particles[index].connector[m]+=(strands-1)*particlePerStrand;
  }
  // box+=minSize;
  return true;
};

bool Analysis::addFalseType(vector<int> colors){
  if(colors.size()==particlePerStrand-1){// As the central particle is always 100 so prevent redundency
    colors.push_back(100);
  } else if(colors.size()<particlePerStrand){// if less number of colors are provided the process will fails
    cout << "Invalid number of colors, operation failed"<<endl;
    return false;
  }
  strands+=1;
  particleTypes+=1;
  particleNum+=particlePerStrand;
  particles.resize(particleNum);
  for(i=0;i<3;i++) box[i]=box[i]*particleTypes/(particleTypes-1); //go over x,y,z of the box
  LR_vector minSize=(box/particleTypes);
  // cout<<minSize.multiplyEach(shift)<<endl;
  for(i=0;i<particlePerStrand;i++){
    int index =(strands-1)*particlePerStrand+i;
    particles[index]=particles[i];
    particles[index].color=colors[i];
  }
  // box+=minSize;
  return true;
};

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

bool Analysis::writePHBTopology(string topology){
  if(topology == "") topology = output + ".top";
  ofstream outputTop(topology, ios::trunc);
  if(!outputTop.is_open()) return false;

  outputTop.precision(15);
  if(patchy) strands=particleNum;
  outputTop<<particleNum<<" "<<strands<<" "<<particleNum<<std::endl;
  for(i=0;i<sourcePatch.size();i++){
    outputTop<<"iP "<<i<<" "<< sourcePatch[i].color<<" "<<sourcePatch[i].strength<<" "<<sourcePatch[i].position.x<<" "<<sourcePatch[i].position.y<<" "<<sourcePatch[i].position.z<<" "<< sourcePatch[i].a1.x<<" "<<sourcePatch[i].a1.y<<" "<<sourcePatch[i].a1.z<<" "<<sourcePatch[i].a2.x<<" "<<sourcePatch[i].a2.y<<" "<<sourcePatch[i].a2.z<<std::endl;
  }
  for(i=0;i<patchConfig.size();i++){
    outputTop<<"iC "<<i<<" ";
    for(int j=0;j<patchConfig[i].size();j++){
      outputTop<<patchConfig[i][j]<<" ";
    }
    outputTop<<std::endl;
  }
  if(patchy){
    for(i=0;i<particleNum;i++){
      outputTop<<"-3 "<<i<<" "<<particles[i].strand<<" "<<patchyRadius;
      outputTop<<std::endl;
    }
  }else{
    for(i=0;i<particleNum;i++){
      outputTop<<-3<<" "<<particles[i].strand<<" "<<particles[i].color<<" "<<particles[i].radius;
      outputTop<<std::endl;
    }
  }
  return true;
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

bool Analysis::populateSingle(int num,double seperator){
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
            particles[index] =particles[p]; //copy the whole particle property
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

bool Analysis::selectIDs(Analysis* selectedParticles,std::vector<int> ids,bool reboxing){
  selectedParticles->particleNum=ids.size();
  selectedParticles->particles.resize(ids.size());
  selectedParticles->box=this->box;
  for(i=0;i<ids.size();i++){
    selectedParticles->particles[i]=this->particles[ids[i]];
  }
  if(reboxing) return selectedParticles->reboxing();
  return true;
}

bool Analysis::reboxing(int offset){
  LR_vector mean = CenterForIndex(-1);
  for(i=0;i<particleNum;i++){
    particles[i].r-=mean;
  }
  //find minimum
  LR_vector min={0,0,0};
  for(i=0;i<particleNum;i++){
    for(int j=0;j<3;j++)
      if(particles[i].r[j]<min[j]) min[j]=particles[i].r[j];
  }
  min-=offset;
  // Adjust the quantity
  for(i=0;i<particleNum;i++){
    particles[i].r-=min;
  }
  inboxing();
  box={0,0,0};
  for(i=0;i<particleNum;i++){
    for(int j=0;j<3;j++)
      if(particles[i].r[j]>box[j]) box[j]=particles[i].r[j];
  }
  for(int j=0;j<3;j++)box[j] = ceil(box[j]);
  return true;
}

bool Analysis::populate(int num, double seperator){
  if(particleTypes==1) return populateSingle(num,seperator);
  if(particlePerStrand==0) particlePerStrand=particleNum/particleTypes;
  if(num%particleTypes!=0) return false;
  cout<<"tell me something"<<endl;

  Analysis psp("","","new"+type); // generating empty particle class
  psp.particles.resize(particleNum);
  for(i=0;i<particleNum;i++){
    psp.particles[i]=particles[i];
  }

  particleNum=particlePerStrand*num;
  strands=num;
  particles.clear(); // flush out everything
  particles.resize(particleNum);

  int dim=std::ceil(std::pow(num,1.0/3.0));
  LR_vector minSize = box-2+seperator;
  // minSize*=1.05;
  box=minSize*dim+2;
  // cout<<psp.particles[0].connector<<endl;
  // for(i=0;i<particleTypes;i++){ //need to work on later to make certain the sample particle is at mindist config
  //   LR_vector minimum = {0,0,0};
  //   for(int j=0;j<particlePerStrand;j++){
  //     if (minimum.x<)
  //   }
  // }
  std::vector<int> count(particleTypes,0); // Keep count of different types of particle, prevent excess
  LR_vector translate={0,0,0};
  int currentNum=1;
  auto drawn=draw();
  // cout<<drawn<<endl;
  buildingParticles:
    for(i=0;i<dim;i++){
      for(int j=0;j<dim;j++){
        for(int k=0;k<dim;k++){
          LR_vector shift={minSize.x*k,minSize.y*j,minSize.z*i};
          for(int p=0;p<particlePerStrand;p++){
            int index=p+k*particlePerStrand+j*dim*particlePerStrand+i*dim*dim*particlePerStrand;
            cout<<index<<" ";
            int pspIndex = p+drawn[currentNum]*particlePerStrand;
            particles[index]=psp.particles[pspIndex];
            particles[index].r+=shift;
            particles[index].strand=(type!="crystal")? currentNum:drawn[currentNum];
            if(type=="crystal") particles[index].patches=psp.particles[pspIndex].patches;
            for(int m=0;m<particles[index].connector.size();m++) particles[index].connector[m]+=(currentNum-1)*particlePerStrand;
          }
          currentNum+=1;
          if(currentNum==num) goto exit;
        }
      }
    }
    return true;
  exit:
    return true;

  cout<<"Program should not have reached this state. Building particles failed."<<endl;
  return false;
}

bool Analysis::PHBhelixExtender(int numParticles,double particleRadius,int numHelix,double helixRadius,LR_vector com){
  PatchesToSpherical(); //converts patches coordinates from cartesian to spherical#
  // Checks for this formulation to work
  if(helixRadius*numHelix>particleRadius){cout<< "Too many big helixes for this algo"<<endl; return false;}
  double xOffset = particleRadius-helixRadius*numHelix; //
  double yGap = (particleRadius*2-helixRadius*2*patchConfig[0].size())/patchConfig[0].size();

  //Set particle parameters
  strands=numParticles;
  particleNum+=(patchConfig[0].size()*numHelix+1)*numParticles;
  particles.resize(particleNum);
  double sourcePatchNum = sourcePatch.size();
  double patchConfigNum = patchConfig.size();
  if (particleTypes != patchConfigNum){
    cout<< "Current Particle Types =" << particleTypes<< "  resetting this"<<endl;
    particleTypes = patchConfigNum;
  }
  if(particlePerStrand==0) particlePerStrand = patchConfig[0].size()*numHelix+1;

 // Setting up Patches configuration
  #pragma omp parallel for
  for(i=0;i<patchConfig[0].size();i++)



  //Real stuff starts here
  #pragma omp parallel for collapse(2)
  for(int j=0;j<patchConfig[0].size();j++){for(i=0;i<numHelix;i++){
    int index = patchConfig[0].size()*j+i;
    particles[index].r={xOffset+helixRadius*2*i,j*2*helixRadius,helixRadius};
    particles[index].radius=helixRadius;
    if(i==0){
      particles[index].connector.push_back(i+1);
    }else if(i==numHelix-1){
      particles[index].connector.push_back(i-1);
      particles[index].r.y+=yGap;
      particles[index+1].r={particleRadius,particleRadius,particleRadius+4*helixRadius};
      particles[index+1].radius= particleRadius;
    }else{
      particles[index].connector.push_back(i+1);
      particles[index].connector.push_back(i-1);
      particles[index].r.y+=yGap;
    }
  }} // This is for the square loop
  // #pragma omp parallel for collapse(2)


  return true;
};



bool Analysis::PatchesToSpherical(){
  #pragma omp parallel for
  for(i=0;i<sourcePatch.size();i++){
    LR_vector converted;
    if(!sourcePatch[i].spherical){
      sourcePatch[i].position = cartesianToSpherical(sourcePatch[i].position);
      sourcePatch[i].spherical=true;
    }
    // #pragma omp critical
    // {
    // cout<<i<<"\t";
    // cout<<sourcePatch[i].position<<endl;
    // }
  }
  return true;
}

vector<int> Analysis::draw(){
  int maxCount=(particleNum/particlePerStrand)/particleTypes;
  vector<int> count(particleTypes,0);
  vector<int> result(strands);
  std::random_device rd;  // a seed source for the random number engine
  std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> distrib(0,particleTypes-1);
  int random;
  for(i=0;i<strands;i++){
    bool cont=true;
    while (cont){
      random = distrib(gen);
      // cout<<random<<endl;
      if(count[random]<maxCount) cont=false;
    }
    count[random]+=1;
    result[i]=random;
  }
  // cout<<"Particle Numbers"<<strands<<endl;
  // cout<<"Count of particle types = "<<count<<endl;
  // cout<<"Length of result = "<<result.size()<<endl;
  return result;
}

std::vector<std::string> npSplit(std::string s, std::string delimiter){
  std::vector<std::string> output;
  size_t pos;
  std::string token;
  while((pos=s.find(delimiter))!=std::string::npos){
    token=s.substr(0,pos);
    istringstream iss(token);
    iss>>token;
    output.push_back(token);
    s.erase(0,pos+delimiter.length());
  }
  output.push_back(s);
  return output;
}

std::vector<std::vector<double> > readCSV(std::string inputFile){
  std::vector<std::vector<double> > output;
  std::string line;
  ifstream input(inputFile);
  while(getline(input,line)){
    std::vector<double> temp;
    auto split = npSplit(line,",");
    for(int i=0;i<split.size();i++)
      temp.push_back(std::stod(split[i]));
    output.push_back(temp);
  }
  return output;
}

bool writeCSV(std::vector<std::vector<double> > data,std::string output="output.csv",int precession=10){
  ofstream outputFile(output);
  outputFile.precision(precession);
  for(int i=0;i<data.size();i++){
    for(int j=0;j<data[i].size();j++){
      if(j==data[i].size()-1){
        outputFile<<data[i][j]<<"\n";
      }else{
        outputFile<<data[i][j]<<",";
      }
    }
  }
  outputFile.close();
  return true;
}

LR_vector cartesianToSpherical(LR_vector cartesian){
  LR_vector spherical;
  spherical.x=cartesian.module();
  spherical.y=acos(cartesian.z/spherical.x);
  spherical.z=atan2(cartesian.y,cartesian.x);
  return spherical;
}

LR_vector sphericalToCartesian(LR_vector spherical){
  LR_vector cartesian;
  cartesian.x=spherical.x*sin(spherical.y)*cos(spherical.z);
  cartesian.y=spherical.x*sin(spherical.y)*sin(spherical.z);
  cartesian.z=spherical.x*cos(spherical.y);
  return cartesian;
}

std::vector<LR_vector> generateRandomColors(int n) {
    if(n<1){std::cout<<"Generate Random Colors says: avoid giving me meaningless input." ;return {};};
    std::vector<LR_vector> colors;
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator
    std::uniform_int_distribution<> distrib(0, 255); // Define the range

    for (int i = 0; i < n; ++i) {
        LR_vector color;
        color.x = distrib(gen); // Generate a random value for red
        color.y = distrib(gen); // Generate a random value for green
        color.z = distrib(gen); // Generate a random value for blue
        colors.push_back(color);
    }

    return colors;
}

std::vector<LR_vector> generateComplimentaryColors(std::vector<LR_vector> colors,double percentage){
  std::vector<LR_vector> output;
  percentage = std::min(percentage,1.0);
  output.resize(colors.size());
  for(int i=0;i<colors.size();i++){
    output[i].x = std::min(255,static_cast<int>(colors[i].x+(255-colors[i].x)*percentage));
    output[i].y = std::min(255,static_cast<int>(colors[i].y+(255-colors[i].y)*percentage));
    output[i].z = std::min(255,static_cast<int>(colors[i].z+(255-colors[i].z)*percentage));
  }
  return output;
}

void skipLines(std::ifstream& file, unsigned int linesToSkip) {
    std::string unused;
    while (linesToSkip-- > 0 && std::getline(file, unused)) {
    }
}