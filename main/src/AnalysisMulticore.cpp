#include "AnalysisMulticore.h"

bool Amc::meanConfig(){
    cout<<"Mean config is called"<<endl;
    // traj.resize(readSegments);
    return true;
}

bool Amc::readPartialTrajectory(ifstream *trajectoryFile,int numRead, int skip){
    if(skip!=0) for(i=0;i<skip*(particleNum+3);i++) getline(*trajectoryFile,line); // Need to check how to jump lines quickly
    for(i=0;i<numRead;i++){
        Traj tempTraj;
        tempTraj.updateParticleNumber(particleNum);
        getline(*trajectoryFile,line);
        tempTraj.time=std::stoi(line);
        getline(*trajectoryFile,line);
        if(i==0){ // Only for the 1st frame set the box, otherwise ignore it.
            ss.clear();ss.str(line);
            ss>>temp; ss>>temp;
            ss>>box.x;ss>>box.y;ss>>box.z;
        }
        getline(*trajectoryFile,line); //Skip Energy
        for(int j=0;j<particleNum;j++){
            ss.clear();ss.str(line);
            ss>>tempTraj.r[j].x;ss>>tempTraj.r[j].y;ss>>tempTraj.r[j].z;
            ss>>tempTraj.a1[j].x;ss>>tempTraj.a1[j].y;ss>>tempTraj.a1[j].z;
            ss>>tempTraj.a3[j].x;ss>>tempTraj.a3[j].y;ss>>tempTraj.a3[j].z;
        }
        traj.push_back(tempTraj);
    }
    return true;
}


// bool Amc::readTrajectory(std::string config,int start, int numRead){
//     ifstream trajectory(config);
//     if(!trajectory.is_open()) return false;
//     int skipline=(particleNum+3)*start;
//     for(i=0;i<skipline;i++) getline(trajectory,line); // This is bad algorithm
//     for(i=0;i<num)
//     return true;
// }

bool Amc::writePHBtopology(std::string topology){
    if (topology == "")
      topology = output + ".top";
    
}

//   bool Analysis::readConfig(string config){
//     ifstream inputConfig(config);
//     if (!inputConfig.is_open())
//       return false;
//     getline(inputConfig, line);
//     getline(inputConfig, line);
//     ss.clear();
//     ss.str(line);
//     ss >> temp;
//     ss >> temp;
//     ss >> box.x;
//     ss >> box.y;
//     ss >> box.z;
//     getline(inputConfig, line);
//     ss.clear();
//     ss.str(line);
//     ss >> temp;
//     ss >> temp;
//     ss >> energy.x;
//     ss >> energy.y;
//     ss >> energy.z;

//     int count=0;
//     if(particleNum==0){
//       particles.resize(99999);
//       while(getline(inputConfig, line)){
//         if(line.empty()||line[0]=='#') continue;
//         ss.clear();
//         ss.str(line);
//         ss >> particles[count].r.x;
//         ss >> particles[count].r.y;
//         ss >> particles[count].r.z;
//         ss >> particles[count].a1.x;
//         ss >> particles[count].a1.y;
//         ss >> particles[count].a1.z;
//         ss >> particles[count].a3.x;
//         ss >> particles[count].a3.y;
//         ss >> particles[count].a3.z;
//         count++;
//       }
//       particles.resize(count);
//       particleNum=count;
//       return true;
//     }

//     for (i = 0; i < particleNum; i++)
//     {
//       getline(inputConfig, line);
//       // if(line.empty()||line[0]=='#') continue;
//       ss.clear();
//       ss.str(line);
//       ss >> particles[i].r.x;
//       ss >> particles[i].r.y;
//       ss >> particles[i].r.z;
//       ss >> particles[i].a1.x;
//       ss >> particles[i].a1.y;
//       ss >> particles[i].a1.z;
//       ss >> particles[i].a3.x;
//       ss >> particles[i].a3.y;
//       ss >> particles[i].a3.z;
//     }
//     return true;
//   }