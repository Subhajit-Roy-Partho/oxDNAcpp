#include "SimulationManager.h"
#include "main.h"
#include "yaml-cpp/yaml.h" //This is an external library
#include <filesystem>
#include <chrono>
#include <thread>
// #include "gnuplot-iostream.h"
using namespace std;
using namespace std::chrono_literals;
    Manager::Manager(string configPath,bool sbatch){
        this->configPath=configPath;
        this->sbatch=sbatch;
        // configuration();
        if(configPath !="") readYAML();
    };

    // bool configuration(){
    //     YAML::Node config = YAML::LoadFile(configPath+"/config.yaml");
    //     projectName = config["project"].as<string>();
    //     outputFiles = config["outputFiles"].as<std::string>();
    //     cout<<projectName<<"\t"<<outputFiles<<"\n";
    //     return true;
    // }

    bool Manager::readYAML(){
        YAML::Node config = YAML::LoadFile(configPath+"/config.yaml");
        replicas = config["replicas"].as<int>();
        variable = config["variable"].as<string>();
        outputFiles = config["outputFiles"].as<string>();
        sbatchLine = config["sbatchLine"].as<int>();
        sbatch = config["sbatch"].as<bool>();
        if(config["exec"]){
            exec=config["exec"].as<string>();
        }
        cout << "This is being executed   "<<exec << endl;

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
        path /= configPath;
        filesystem::current_path(path);
        path = filesystem::current_path();
        return true;
    };

    bool Manager::setup(){
        // auto path = filesystem::current_path();
        // path /= configPath;
        // filesystem::current_path(path);
        // path = filesystem::current_path();

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
                    if(sbatch){
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
                        tempInt=system("sbatch submit.sh");
                    }else{
                        filesystem::current_path(relative);
                        tempString=exec+" input >> out.txt &";
                        tempInt=system(tempString.c_str());
                    }

                    allocations[currentCluster]-=1;
                }
            }
        }

        return true;
    };



    bool Manager::plot(){
        // std::cout <<"this is working"<<std::endl;
        cout<<path<<endl;
        std::ofstream commandfile;
        for (int p=0;p<inputFiles.size();p++){
            for(int j=0;j<replicas;j++){
                commandfile.open("plot.gnuplot",ofstream::out|ofstream::trunc);
                    commandfile<<"set term png size 1200,900\n";
                    commandfile<<"set xlabel 'Time(SU)'\n";
                    commandfile<<"set logscale x\n";
                    commandfile<<"set yrange [0:-1.4]\n";
                    commandfile<<"set ylabel 'Energy(SU)'\n";
                    commandfile<<"set title '"+projectName[p]+"'\n";
                    commandfile<<"set output './plots/"+projectName[p]+"-r"+to_string(j+1)+".png'\n";
                for (int t=0;t<values.size();t++){
                    string relative=outputFiles+"/"+projectName[p]+"/replica-"+to_string(j+1)+"/"+string(variable)+"-"+to_string_with_precision(values[t],3)+"/energy.dat";
                    // cout << relative<<endl;
                    if(t==0){
                        commandfile<<"plot '"+relative+"' w l title '"+string(variable)+"-"+to_string_with_precision(values[t],3)+"',\\\n";
                    }else if(t==values.size()-1){
                        commandfile<<"\t'"+relative+"' w l title '"+string(variable)+"-"+to_string_with_precision(values[t],3)+"'\n";
                    }else{
                        commandfile<<"\t'"+relative+"' w l title '"+string(variable)+"-"+to_string_with_precision(values[t],3)+"',\\\n";
                    }
                }
                commandfile.close();
                tempInt=system("gnuplot -p plot.gnuplot");
                // this_thread::sleep_for(2s);
            }
        }
        cout<<"Done plotting\n";
        removeList.push_back("plot.gnuplot");
        return true;
    }

    // bool clean(){
        
    // }