// To generate crystal seeds
// crystal outputName topfile
//Error 10 = can't open top file

#include "def.cpp"
using namespace std;

int main(int argNum, char *argv[]){
    
    Analysis crystal("BIG.top","BIG.conf","crystal");
    int clusters[] = {1070,134,145,155,1004,1385,560,1291,904,136,1112,126,1355,1850,1384,686,392,1605,618,1024,2027,474,1640,1720,1831,387,1504,947,1832,1262,676,1067,195,121,734,245};
    
    crystal.pickAndPlace(clusters);
    
    return 0;

}