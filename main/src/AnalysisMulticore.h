#include "Analysis.h"

class Amc:Analysis{
public:
    int readSegments=20;
    bool meanConfig();
    bool readTrajectory(std::string config="",int start=0, int numRead=1);
    bool readPartialTrajectory(ifstream *trajectoryFile,int numRead,int skip=0);
    bool writePHBtopology(std::string topology="");
};