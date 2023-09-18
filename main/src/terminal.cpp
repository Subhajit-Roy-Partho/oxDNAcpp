/*
    Usage: terminal manager --config "location"
*/

#include "Analysis.cpp"
#include "SimulationManager.cpp"


int main(int argc, char *argv[]){
    // Analysis randPart("","","newCrystal");



    if(!std::strcmp(argv[1],"manager")){
        for(int i=2;i<argc;i++){
            if(!std::strcmp(argv[i],"--config")){
                Manager manage(std::string(argv[i+1]));
                manage.readYAML();
            };
        }
    }


    // std::cout<< "working"<<std::endl;
    return 0;
}