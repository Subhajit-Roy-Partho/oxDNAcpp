/* ############## Errors #####################
10 - error opening input Topology
11 - error opening input Config
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "LR_vector.h"
#include <gsl/gsl_linalg.h>

// int gsl_linalg_SV_decomp(gsl_matrix *A, gsl_matrix *V, gsl_vector *S, gsl_vector *work)

using namespace std;

class Particle{
    public:
        int id,color,strand;
        LR_vector r,a1,a3;
};

class Analysis{
    public:
        int particleNum,strands,i;
        string type;
        LR_vector box,energy;
        std::vector<Particle> particles;
        gsl_matrix *V=gsl_matrix_alloc(3,3);
        gsl_vector *S = gsl_vector_alloc(3);
        gsl_vector *work = gsl_vector_alloc(3);

        Analysis(string topology,string config,string type="", string output="output",string externalForces="",string parameter1="",string parameter2=""){
            if(type=="crystal"){
                this->type=type;
                readCrystalTopology(topology);
                readConfig(config);
                // pickAndPlace();
            }
        }
        ~Analysis(){
            gsl_matrix_free(V);
            gsl_vector_free(S);
            gsl_vector_free(work);
        }

        bool pickAndPlace(int *cluster){
                

            double points[] = {158720.15575206,42724.03921793,56622.47200362,
                            42724.03921793,132381.4182789,-83288.45034046,
                            56622.47200362,-83288.45034046,100855.62941231};
            gsl_matrix_view U = gsl_matrix_view_array(points,3,3);
            int pass = gsl_linalg_SV_decomp(&U.matrix,V,S,work);
            gsl_matrix_fprintf(stdout,V,"%g");
            return true;
        };

        bool planeFitting(int *cluster, LR_vector center, LR_vector normal){
            
            return true;
        }

    private:
        string line,temp;
        istringstream ss;
        int readCrystalTopology(string topology){
            ifstream inputTop(topology);
            if(!inputTop.is_open()) return 10;
            getline(inputTop,line);
            ss.clear();
            ss.str(line);
            ss>>particleNum;
            ss>>strands;
            particles.resize(particleNum);

            ss.clear();
            getline(inputTop,line);
            ss.str(line);

             for(i=0;i<particleNum;i++){
                ss>>temp;
                particles[i].id=i;
                particles[i].color=stoi(temp);
            }
            return 0;
        }

        int readConfig(string config){
            ifstream inputConfig(config);
            if(!inputConfig.is_open()) return 11;
            getline(inputConfig,line);
            getline(inputConfig,line);
            ss.clear();ss.str(line);ss>>temp; ss>>temp;
            ss>>box.x; ss>> box.y;ss>>box.z;
            getline(inputConfig,line);
            ss.clear();ss.str(line);ss>>temp;ss>>temp;
            ss>>energy.x;ss>>energy.y;ss>>energy.z;
            for(i=0;i<particleNum;i++){
                getline(inputConfig,line);
                ss.clear();
                ss.str(line);
                ss>>particles[i].r.x;ss>>particles[i].r.y;ss>>particles[i].r.z;
                ss>>particles[i].a1.x;ss>>particles[i].a1.y;ss>>particles[i].a1.z;
                ss>>particles[i].a3.x;ss>>particles[i].a3.y;ss>>particles[i].a3.z;
            }
            return 0;
        }

        int readPatches(string patches){
            ifstream inputPatches(patches);
            return 0;
        }

        int readParticles(string crystalpar){
            ifstream inputCrystal(crystalpar);
            return 0;
        }
};

// class np{
//     public:
//     <template typename A>
//     void mean(){

//     }
// };

// Selected Points
// 1070,134,145,155,1004,1385,560,1291,904,136,1112,126,1355,1850,1384,686,392,1605,618,1024,2027,474,1640,1720,1831,387,1504,947,1832,1262,676,1067,195,121,734,245,217