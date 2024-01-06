#include "Analysis.h"
#include "omp.h"

class Amc:Analysis{
public:
    int readSegments=20;
    bool meanConfig();
    bool readTrajectory(std::string config="",int start=0, int numRead=1);
    bool readPartialTrajectory(ifstream *trajectoryFile,int numRead,int skip=0);
};