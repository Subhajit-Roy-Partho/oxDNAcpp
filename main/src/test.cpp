#include <iostream>
// #include "Analysis.h"
#include "AnalysisMulticore.h"
using namespace std;

// Return selected patchy particles
void patchyReturn(){
    Analysis crystals("../../managerExample/4flat/input.top","../../managerExample/4flat/input.dat","crystal");
    std::vector<int> ids ={1174,1168,1942,620};
    Analysis selected;
    crystals.selectIDs(&selected,ids);
    selected.writeCrystalTopology();
    selected.writeConfig();
};

void MGLgenerator(){
    //Analysis crystals("../../managerExample/4flat/input.top","../../managerExample/4flat/input.dat","crystal");
    Analysis crystals("../../managerExample/4flat/input.top","/tmp/subho/last_conf-t125-b1000.dat","crystal");
    crystals.readCrystalPatches("../../managerExample/4flat/sat24.patches.txt");
    crystals.readCrystalParticlePatchyConfig("../../managerExample/4flat/CRYSTAL.particles.txt");
    crystals.writeMGL("output.mgl",0.4,0.1);
    // cout<<crystals.particles[500].patches<<endl;
}

void PSPgenerator(){
    Analysis ico("../../managerExample/DNA/icoAligned.top","../../managerExample/DNA/icoAligned.dat","DNA"); // Loading the top and the dat files
    // ico.shiftbox({0,0,0}); // Bringing everything to 1st coordinate
    // ico.testBoxOverloaded(); // Checking overload from previous action.

    vector<vector<int>> ids{
        {6532,6555,6559,6589,6158,6217,4906,6182,6180,4925}, //Red 100
        {4471,4492,4497,5028,5033,5054,5059,5080,5085,4466}, //Blue - purple-5
        {4023,4044,4050,4070,6013,6033,6359,6380,6406,6386}, //Yellow - turcoise-6
        {2205,6072,6077,6099,5232,5211,2515,6512,5257,2826}, //Green - 100
        {2566,2781,3668,3497,3522,3549,3529,3554,3637,2722}, // Purple - blue-2
        {4678,1378,4755,4705,4781,4762,1484,1529,4672,4814}, // Turcoise - 100
        {4352,4179,4186,4207,4213,4294,4300,4320,4326,2004}, // Pink - green-4
        {5389,746,695,128,5422,5442,5306,5449,5312,5333},//Brown - 100
        {3186,3207,3233,3213,3239,2885,2877,3378,3384,3180}, //white - red -1
        {6600,6620,275,6646,6314,6293,390,436,5735,5714}, // brown - yellow-3
        {2360,2318,2311,2264,2256,2471,2421,2461,3773,3794}, //100
        {6225,6246,224,180,171,652,6665,6687,589,6712} //100
        };
    
    vector<int> colors{100,25,26,100,22,100,24,100,21,23,100,100}; // if 12 colors provided 13th color will be by default 100 or colorless
    vector<double> radius{8,10}; // if 2 of the radius are provided then 1st for all particle and 2nd for the central particle
    // cout << ids.size()<<endl;
    Analysis psp("","","newPSP"); // generating empty particle class
    ico.generatePSP(&psp,ids,colors,radius,5,100); // generating ccg particle from oxDNA file.
    colors ={100,-25,30,100,-22,100,29,100,27,28,100,100};
    psp.addFalseType(colors);
    colors = {100,-30,-24,100,-28,100,32,100,31,-21,100,100};
    psp.addFalseType(colors);
    colors = {100,-32,23,100,-31,100,-29,100,-27,-23,100,100};
    psp.addFalseType(colors);
    psp.populate(128,5); //generate 125 crystal with extra 1 su seperation between them
    // psp.reboxing();
    // psp.boxToCubic();// make the box cubic nothing is effected
    psp.writeCCGtopology("127.top"); // write the CCG topology
    psp.writeCCGviewTopology("127view.top"); // write oxview readeable topology
    psp.writeConfig("127.dat"); // write config dat file
    // cout <<psp.particles[5].connector<<endl;
};

void CubeGenerator(){
    Analysis cub("","","newcrystal");
    cub.particleNum=2,cub.particlePerStrand=1,cub.strands=2; cub.particles.resize(2),cub.particleTypes=2;
    cub.box=(LR_vector){1.2,1.2,1.2};
    Patch cubSide;
    cubSide.position=(LR_vector){0,0.5,0};cubSide.color=22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){0.5,0,0};cubSide.color=22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){0,-0.5,0};cubSide.color=22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){-0.5,0,0};cubSide.color=22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){0,0,-0.5};cubSide.color=22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){0,0,0.5};cubSide.color=22;cub.sourcePatch.push_back(cubSide);

    //complimentary patches
    cubSide.position=(LR_vector){0,0.5,0};cubSide.color=-22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){0.5,0,0};cubSide.color=-22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){0,-0.5,0};cubSide.color=-22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){-0.5,0,0};cubSide.color=-22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){0,0,-0.5};cubSide.color=-22;cub.sourcePatch.push_back(cubSide);
    cubSide.position=(LR_vector){0,0,0.5};cubSide.color=-22;cub.sourcePatch.push_back(cubSide);

    cub.particles[0].patches={0,1,2,3,4,5}; cub.particles[0].strand=0;
    cub.particles[1].patches={6,7,8,9,10,11};cub.particles[1].strand=1;
    cub.populate(1000,3);
    for(int i=0;i<cub.particleNum;i++){
        cout<<cub.particles[i].patches<<endl;
    }
    cout<<cub.type<<endl;
    cub.writeMGL("output.mgl",0.4,0.1);
}

void PatchyGenerator(){
    
}

void DNAarmAnalysis(){
    Analysis dna("../../managerExample/DNA/icoAligned.top","../../managerExample/DNA/icoAligned.dat","DNA");
}

void PHBgenerator(){
    Analysis phb("newPHB");
    phb.readCrystalPatches("../../managerExample/4flat/sat24.patches.txt");
    phb.readCrystalParticlePatchyConfig("../../managerExample/4flat/CRYSTAL.particles.txt");
    // phb.PatchesToSpherical();
    // cout<<sphericalToCartesian(LR_vector{0.5,M_PI_2,M_PI})<<endl;
    phb.PHBhelixExtender();
    // phb.PHBhelixExtender();
}

void patchyToPHB(){
    Analysis crystals("../../managerExample/4flat/input.top","../../managerExample/4flat/input.dat","crystal");
    crystals.readCrystalPatches("../../managerExample/4flat/sat24.patches.txt");
    crystals.readCrystalParticlePatchyConfig("../../managerExample/4flat/CRYSTAL.particles.txt");
    crystals.patchy=true;
    crystals.writePHBTopology("output.top");
    crystals.writeConfig("output.dat");
}

int main(){
    // patchyReturn();
    // PSPgenerator();
    MGLgenerator();
    // Analysis("","","newcrystal");
    // CubeGenerator();
    // PHBgenerator();
    // patchyToPHB();

    return 0;
}

//PHB system 
// 55 ico, 39 helix ratio(helix,ico)=0.71