// To generate crystal seeds
// crystal outputName topfile
//Error 10 = can't open top file

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "def.cpp"
using namespace std;

int main(int argNum, char *argv[]){
    string line,temp;
    int particleNum,types,i=0;

    ifstream inputTop(argv[2]);
    if(!inputTop.is_open()) return 10;

    
    getline(inputTop,line);
    istringstream ss (line);
    ss>>particleNum;
    ss>>types;
    auto particles = new Particle[particleNum];

    getline(inputTop,line);
    ss.clear();
    ss.str(line);

    for(i=0;i<particleNum;i++){
        ss>>temp;
        particles[i].id=i;
        particles[i].color=stoi(temp);
    }


}