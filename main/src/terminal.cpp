/*
    Usage: terminal manager --config "location in folder"
*/

// #include "main.h"
#include "SimulationManager.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[]){
    string config="";
    bool sbatch=true;
    for(int i=2;i<argc;i++){
        if(string(argv[i])=="--config"){
            config=string(argv[i+1]);
        }
        if(string(argv[i])=="--sbatch"){
            sbatch=bool(argv[i+1]);
        }
    }

    Manager manage(config);




    if(string(argv[1])=="manager")manage.setup();
    if(string(argv[1])=="plot")manage.plot();


    return 0;
}