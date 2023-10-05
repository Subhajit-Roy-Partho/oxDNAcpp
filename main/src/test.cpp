#include <iostream>
#include "Analysis.h"
using namespace std;

int main(){
    Analysis ico("../../managerExample/DNA/ico.top","../../managerExample/DNA/ico.dat","DNA");
    cout <<ico.shiftbox({1,0,0})<<endl;
    return 0;
}