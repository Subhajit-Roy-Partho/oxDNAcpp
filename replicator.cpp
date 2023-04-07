/* Error Code:
    5 - topology input file not found.
    6 - error in opening output topology file.
    52 - error opening input coordinate file.
*/

// Usage replicator topologyfile replica distance coordfile

// cg = type, total count, joints

#include <iostream>
#include <fstream>
#include <string>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <vector>
#include <cmath>
#include <iomanip>
// #include <jsoncpp/json/json.h>
// #include<jsoncpp/json/value.h>

using namespace std;
using namespace xt;



int main(int argc, char* argv[]){
    double distance = 70; // argc[3]

// Topological File

    xarray<int> cg = zeros<int>({1000,30});
    string line,name[1000],temp,header; 
    int j,number,type,i=0,count=0,countN=0,lastcount=0,numcount=0,templast=0;
    ifstream topfile (argv[1]);
    if (topfile.is_open()){
        getline(topfile,line);
        header = line;
        while(getline(topfile,line)){
            j=0;type=0;
            istringstream ss (line);
            while(ss.tellg() !=-1){
                ss>>temp;
                if(j==1){
                    name[i] = temp;
                }else{
                    cg(i,j)=stoi(temp);
                }
                if(j==0 && stoi(temp)<0)
                    count++;
                if(j==0 && stoi(temp)>0)
                    countN++;
                j++;
            }
            cg(i,1)=j;
            i++;
        }
    }else{
        cout << "Input file not found"<<endl;
        return 5;
    }
    topfile.close();
    int replicas = stoi(argv[2]);
    // int replicas =27;
    // cout << replicas <<"\n";

    ofstream outfile ("replica.top");
    if(!outfile.is_open()){
        cout<<"Error in opening output file";
        return 10;
    }
    // outfile<<header<<"\n";
    istringstream ss (header);
    while(ss.tellg() !=-1){
        ss>>temp;
        outfile<< stoi(temp)*replicas<<" ";
    }
    outfile<<endl;
    for(int i=0;i<replicas;i++){
        for(int j=0;j<count;j++){
            for(int p=0;p<cg(j,1);p++){
                if(p==1){
                    outfile<<name[j]<<" ";
                }else if(p==0){
                    outfile<<cg(j,0)-i<<" ";
                }else {
                    if(cg(j,p) < 0){
                        outfile<<cg(j,p)<<" ";
                    }else{outfile<<cg(j,p)+count*i<<" ";}
                }
            }
            outfile<<"\n";
        }
    }
    lastcount = (replicas-1)*count;

    for(int i=0;i<replicas;i++){
        for(int j=count;j<countN+count;j++){
            for(int p=0;p<cg(j,1);p++){
                if(p==0){
                    outfile<<cg(j,0)+numcount<<" ";
                }else if(p==1){
                    outfile<<name[j]<<" ";
                }else{
                    if(cg(j,p) < 0){
                        outfile<<cg(j,p)<<" ";
                    }else{
                        outfile<<cg(j,p)+lastcount<<" ";
                        templast = (templast <cg(j,p)+lastcount)?cg(j,p)+lastcount:templast;
                    }
                }
            }
            outfile<<"\n";
        }
        numcount+=cg(countN+count-1,0);
        lastcount = templast-count+1;
    }
    outfile.close();
    // cout << row(cg,0)<<endl;

//////////////////////// Data File Modification  /////////////////////////////////////////////////////////
    cout <<count<<endl;
    xarray<float> cgdat = zeros<float>({1000,15});
    i=0,j=0;
    ifstream datfile ("output.dat");
    if(!datfile.is_open()){
        cout << "Error opening input datafile"<<endl;
        return 52;
    }
    getline(datfile,line);
    string timeHeader= line;
    getline(datfile,line);
    string dimHeader=line;
    getline(datfile,line);
    string eHeader = line;
    while(getline(datfile,line)){
        istringstream ss (line);
        for(int dat=0;dat<15;dat++){
            ss>>temp;
            cgdat(i,dat)=stof(temp);
        }
        i++;
    }
    datfile.close();

    ofstream odat ("replica.dat");
    if(!odat.is_open()){
        cout << "Error creating output datafile"<<endl;
        return 62;
    }
    int topP=ceil(pow(replicas, 1/3.));
    int perm[topP*topP*topP][3];
    int pcount=0;
    for(int i=0;i<topP;i++){
        for(int j=0;j<topP;j++){
            for(int p=0;p<topP;p++){
                perm[pcount][0]=i;
                perm[pcount][1]=j;
                perm[pcount][2]=p;
                pcount++;
            }
        }
    }
    odat<<timeHeader<<"\n";
    // odat<<dimHeader<<"\n";
    odat<<"b = 500 500 500\n";
    odat<<eHeader<<"\n";
    float move[3] = {0,0,0};
    float space = 70.0;
    for(int rep=0;rep<replicas;rep++){
        for(int i=0;i<count;i++){
            for(int j=0;j<15;j++){
                if (j<3){ odat<<cgdat(i,j)+perm[rep][j]*space<<" ";}else{odat<<cgdat(i,j)<<" ";}
            }
            odat<<"\n";
        }
        
    }
    move[0]=0;
    for(int rep=0;rep<replicas;rep++){
        for(int i=count;i<count+countN;i++){
            for(int j=0;j<15;j++){
                if (j<3){ odat<<cgdat(i,j)+perm[rep][j]*space<<" ";}else{odat<<cgdat(i,j)<<" ";}
            }
            odat<<"\n";
        }
        move[0]=100;
    }
    odat.close();
    cout<< count<<endl;

//////////////////////////////////////// Parameter File //////////////////////////////////////////////////////////
    xarray<double> cgpar = zeros<double>({1000,5});
    i=0;j=0;
    ifstream par ("output.par");
    if(!par.is_open()){
        cout << "Error opening input Parameter File"<<endl;
        return 53;
    }
    getline(par,line);
    while(getline(par,line)){
        istringstream ss (line);
        for(int dat=0;dat<5;dat++){
            ss>>temp;
            if(dat==3){
                cgpar(i,dat)= (temp =="s") ? 1 : 2;
            }else{
                cgpar(i,dat)=stod(temp);
            }
        }
        i++;
    }
    par.close();


    ofstream opar ("replica.par");
    j=0;
    if(!opar.is_open()){
        cout << "Error creating output parameter file"<<endl;
        return 63;
    }
    opar<<count*replicas<<"\n";
    for(int rep=0;rep<replicas;rep++){
        for(int p=0;p<i;p++){
            for(int j=0;j<5;j++){
                if(j==3){
                    char s = (cgpar(p,j)==1)?'s':'o';
                    opar<<setprecision(17)<<s<<" ";
                }else if(j==0 || j==1){
                    opar<<cgpar(p,j)+rep*count<<" ";
                }else{
                    opar<<cgpar(p,j)<<" ";
                }
            }
            opar<<"\n";
        }
    }
    opar.close();
    cout <<count<<endl;
/////////////////////////////////////////// Force File ////////////////////////////////////////////////////
    int particles[1000],refParticles[1000],realParticle,realRefParticle;
	ifstream force ("external_forces.txt");
	if(!force.is_open()){
        	cout << "Error opening input Force File"<<endl;
        	return 53;
   	}
    // Json::Value jsonvalue;
    // Json::Reader reader;
    // reader.parse(force,jsonvalue);
    // cout << "Total json data"<<jsonvalue <<endl;
    
    i=0;
    while (getline(force,line)){
        istringstream ss (line);
        ss>>temp;
        // cout << temp <<endl;
        if(temp.compare("particle")==0){
            ss>>temp;
            ss>>temp;
            particles[i] = stoi(temp);
        }
        if(temp.compare("ref_particle")==0){
            ss>>temp;
            ss>>temp;
            refParticles[i]=stoi(temp);
            i++;
        }
    }
    
    force.close();
    ofstream oforce ("replica_forces.txt");
    if(!oforce.is_open()){
        cout << "Error creating output Force File"<<endl;
        return 53;
   	}
    // for(int rep=0;rep<replicas;rep++){
    //     for(int j=0;j<i;j++){
    //         // realParticle = (particles[j]<count)?particles[j]+count*(replicas-1):particles[j]+count*(replicas-1)+(count+countN)*(replicas-1);
    //         // realRefParticle = (refParticles[j]<count)?refParticles[j]+count*(replicas-1):refParticles[j]+count*(replicas-1)+(count+countN)*(replicas-1);
    //         oforce<<"{\n\ttype = mutual_trap\n\tparticle = "<<realParticle<<"\n\tref_particle = "<<realRefParticle<<"\n\tstiff = 2.8\n\tr0 = 1.2\n\tPBC = 1\n}"<<endl;
    //     }
    // }

    for(int rep=0; rep<replicas;rep++){
        for(int j=0;j<i;j++){
            int par = (particles[j]<count)?particles[j]+count*(rep):(count*replicas)+(countN*rep)+(particles[j]-count);
            int refPar = (refParticles[j]<count)?refParticles[j]+count*(rep):(count*replicas)+(countN*rep)+(refParticles[j]-count);
            oforce<<"{\n\ttype = mutual_trap\n\tparticle = "<<par<<"\n\tref_particle = "<<refPar<<"\n\tstiff = 2.8\n\tr0 = 1.2\n\tPBC = 1\n}"<<endl;
        }
    }
    oforce.close();
    
//Return no error
    return 0;
}
