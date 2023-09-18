#include <iostream>
#include "yaml-cpp/yaml.h" //This is an external library
using namespace std;
class Manager{
public:
    string configPath,projectName;
    Manager(string configPath=""){
        this->configPath=configPath;
        if(configPath !="") readYAML();
    };

    bool readYAML(){
        YAML::Node config = YAML::LoadFile(configPath);
        projectName = config["project"].as<string>();
        cout<<projectName<<"\n";
        return true;
    };
};