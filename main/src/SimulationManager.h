#include "yaml-cpp/yaml.h" //This is an external library
#include <filesystem>
#include <chrono>
#include <thread>
// #include "gnuplot-iostream.h"
using namespace std;
using namespace std::chrono_literals;
class Manager{
public:
    std::string configPath,variable,outputFiles,line;
    int currentCluster=0,sbatchLine,replicas;
    vector<string> account,queue,inputFiles,projectName,removeList;
    vector<int> allocations;
    vector<double> values;
    YAML::Node temp;
    filesystem::path path = filesystem::current_path();
    bool sbatch;
    Manager(string configPath="",bool sbatch=true);
    bool readYAML();
    bool setup();
    bool plot();
};