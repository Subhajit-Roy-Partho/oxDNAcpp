/* Error Code:
    5 - file not found
*/

// Usage replicator topologyfile

#include <iostream>
#include <fstream>
#include <string>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <vector>
#include <json/json.h>
using namespace std;
using namespace xt;

void PermGenerator(int n, int k)
{
    std::vector<int> d(n);
    std::iota(d.begin(),d.end(),1);
    cout << "These are the Possible Permutations: " << endl;
    do
    {
        for (int i = 0; i < k; i++)
        {
            cout << d[i] << " ";
        }
        cout << endl;
        std::reverse(d.begin()+k,d.end());
    } while (next_permutation(d.begin(),d.end()));
}

int main(int argc, char* argv[]){
    int number[] = {0,0,0};
    int count=0;
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            for(int p=0;p<3;p++){
                cout<< i<<j<<p<<endl;
                count++;
            }
        }
    }
    cout<<count<<endl;
}
