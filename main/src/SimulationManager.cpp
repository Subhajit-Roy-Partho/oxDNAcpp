#include <iostream>
#include <vector>
#include "yaml-cpp/yaml.h" //This is an external library
#include <filesystem>
#include "main.h"
using namespace std;
class Manager{
public:
    std::string configPath,variable,outputFiles,line;
    int currentCluster=0,sbatchLine,replicas;
    vector<string> account,queue,inputFiles,projectName;
    vector<int> allocations;
    vector<double> values;
    YAML::Node temp;
    bool sbatch;
    Manager(string configPath="",bool sbatch=true){
        this->configPath=configPath;
        this->sbatch=sbatch;
        // configuration();
        if(configPath !="") if(readYAML()) setup();
    };

    // bool configuration(){
    //     YAML::Node config = YAML::LoadFile(configPath+"/config.yaml");
    //     projectName = config["project"].as<string>();
    //     outputFiles = config["outputFiles"].as<std::string>();
    //     cout<<projectName<<"\t"<<outputFiles<<"\n";
    //     return true;
    // }

    bool readYAML(){
        YAML::Node config = YAML::LoadFile(configPath+"/config.yaml");
        replicas = config["replicas"].as<int>();
        variable = config["variable"].as<string>();
        outputFiles = config["outputFiles"].as<string>();
        sbatchLine = config["sbatchLine"].as<int>();

        cout<<projectName<<"\n";

        temp = config["project"];
        for(size_t i=0; i<temp.size();i++){
            projectName.push_back(temp[i].as<string>());
        }

        temp = config["inputFiles"];
        for(size_t i=0; i<temp.size();i++){
            inputFiles.push_back(temp[i].as<string>());
        }

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

        for(int p=0;p<projectName.size();p++){
            for(int j=0;j<replicas;j++){
                filesystem::current_path(path);
                string relative2 = "replica-"+to_string(j+1);
                relative2 = outputFiles+"/"+projectName[p]+"/"+relative2;
                for(int i=0;i<values.size();i++){
                    if(allocations[currentCluster]==0){
                        currentCluster+=1;
                        if (currentCluster>account.size()){
                            cout<< "exhausted all the clusters";
                            exit(5);
                        }
                    }
                    filesystem::current_path(path);
                    string relative = string(variable)+"-"+to_string_with_precision(values[i],3);
                    relative=relative2+"/"+relative;
                    filesystem::create_directories(relative);
                    filesystem::copy(inputFiles[p],relative);
                    ofstream file(relative+"/input",ios::out|ios::app);
                    file<< variable<<" = "<<values[i];
                    file.close();
                    file.open(relative+"/submit.sh");
                    ifstream file2(inputFiles[p]+"/submit.sh");
                    for(int i=0;i<sbatchLine-1;i++){
                        getline(file2,line);
                        file<<line<<"\n";
                    }
                    file<<"#SBATCH -p "<<account[currentCluster]<<"\n";
                    file<<"#SBATCH -q "<<queue[currentCluster]<<"\n";
                    while(getline(file2,line)){
                        file<<line<<"\n";
                    }
                    file.close();
                    file2.close();
                    filesystem::current_path(relative);
                    // system("sbatch submit.sh");
                    // if(sbatch){
                    //     ofstream 
                    // }

                    allocations[currentCluster]-=1;
                }
            }
        }

        return true;
    };
};