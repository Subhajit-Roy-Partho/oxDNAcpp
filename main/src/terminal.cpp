/*
    Usage: terminal manager --config "location in folder"
*/

// #include "Analysis.cpp"
#include "SimulationManager.cpp"
#include <iostream>
using namespace std;

int main(int argc, char *argv[]){
    // Analysis randPart("","","newCrystal");
    // cout <<argv[2];



    if(string(argv[1])=="manager"){
        for(int i=2;i<argc;i++){
            if(string(argv[i])=="--config"){
                Manager manage(string(argv[i+1]));
            };
        }
    }


    // std::cout<< "working"<<std::endl;
    return 0;
}