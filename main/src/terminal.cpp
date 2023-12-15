/*
    Usage: terminal manager --config "location in folder"
*/

// #include "main.h"
#include "SimulationManager.h"
#include "Analysis.h"
#include <iostream>
#include <filesystem>
#include <fstream>
using namespace std;
namespace fs = std::filesystem;

std::string extensionFileSearch(std::string path,std::string extension){
    for(auto &p:fs::recursive_directory_iterator(path))
        if(p.path().string().find(extension)!=std::string::npos) 
            return p.path().string();
    return "failed";
};

int main(int argc, char *argv[]){
    string config="",inputTop="input.top",outputTop="output.top",inputDat="input.dat",outputDat="output.dat",correctDat="correct.dat";
    bool sbatch=true;
    for(int i=2;i<argc;i++){
        if(string(argv[i])=="--config"){
            config=string(argv[i+1]);
        }
        if(string(argv[i])=="--sbatch"){
            sbatch=bool(argv[i+1]);
        }
        if(string(argv[i])=="--inputTop"){
            inputTop=string(argv[i+1]);
        }
        if(string(argv[i])==("--inputDat")){
            inputDat=string(argv[i+1]);
        }
        if(string(argv[i])=="--outputDat"){
            outputDat=string(argv[i+1]);
        }
        if(string(argv[i])=="--correctDat"){
            correctDat=string(argv[i+1]);
        }
    };
    if(string(argv[1])=="test"){
        Manager manage(config);
        // manage.readYAML();
    }
    if(string(argv[1])=="inboxDat"){
        Analysis data("",inputDat);
        data.readConfig(inputDat);
        data.inboxing();
        data.writeConfig(outputDat);
    }
    if(string(argv[1])=="manager" || string(argv[1])=="plot"){
        Manager manage(config);
        if(string(argv[1])=="manager")manage.setup();
        if(string(argv[1])=="plot"){cout<<"Getting plots ready for you"<<endl;manage.plot();}
    }else if(string(argv[1])=="analysis"){
    }else if(string(argv[1])=="correctA"){
        Analysis crystal(inputTop,inputDat, "crystal");
        crystal.correctA(correctDat);
        crystal.writeConfig(outputDat);
        return 0;
    }
    else if(std::string(argv[1])=="convert"){
        if(outputDat=="output.dat") outputDat="output.mgl";
        auto topfile = extensionFileSearch(config,".top");
        auto datfile = extensionFileSearch(config,".dat");
        auto particleFile = extensionFileSearch(config,".particles.txt");
        auto patchesFile = extensionFileSearch(config,".patches.txt");

        Analysis crystals(topfile,datfile,"crystal");
        crystals.readCrystalPatches(patchesFile);
        crystals.readCrystalParticlePatchyConfig(particleFile);
        crystals.writeMGL(outputDat,0.4,0.1);
    };


    return 0;
}
