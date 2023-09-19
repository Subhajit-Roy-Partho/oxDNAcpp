#include <iostream>
#include <vector>
#include "yaml-cpp/yaml.h" //This is an external library
#include <filesystem>
#include "main.h"
using namespace std;
class Manager{
public:
    string configPath,projectName,variable;
    vector<string> account,queue;
    vector<int> allocations;
    vector<double> values;
    YAML::Node temp;
    bool sbatch;
    int replicas;
    Manager(string configPath="",bool sbatch=true){
        this->configPath=configPath;
        this->sbatch=sbatch;
        if(configPath !="") if(readYAML()) setup();
    };

    bool readYAML(){
        YAML::Node config = YAML::LoadFile(configPath+"/config.yaml");
        projectName = config["project"].as<string>();
        replicas = config["replicas"].as<int>();
        variable = config["variable"].as<string>();
        cout<<projectName<<"\n";


        temp = config["account"];
        for(size_t i=0; i<temp.size();i++){
            account.push_back(temp[i].as<string>());
        }
        temp = config["queue"];
        for(size_t i=0; i<temp.size();i++){
            queue.push_back(temp[i].as<string>());
        }
        temp=config["allocations"];
        for(size_t i=0; i<temp.size();i++){
            allocations.push_back(temp[i].as<int>());
        }

        temp=config["values"];
        for (size_t i=0; i<temp.size();i++){
            values.push_back(temp[i].as<double>());
        }
        // cout<<values;
        return true;
    };

    bool setup(){
        auto path = filesystem::current_path();
        path /= configPath;
        filesystem::current_path(path);
        path = filesystem::current_path();

        for(int j=0;j<replicas;j++){
            string relative2 = "replica-"+to_string(j+1);
            relative2 = projectName+"/"+relative2;
            for(int i=0;i<values.size();i++){
                string relative = string(variable)+"-"+to_string_with_precision(values[i],3);
                relative=relative2+"/"+relative;
                filesystem::create_directories(relative);
                filesystem::copy("inputs/",relative);
                ofstream file(relative+"/input",ios::out|ios::app);
                file<< variable<<" = "<<values[i];
                file.close();
                // if(sbatch){
                //     ofstream 
                // }
            }
        }

        return true;
    };
};