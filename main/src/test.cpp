#include <iostream>
#include "Analysis.h"
using namespace std;

int main(){
    Analysis ico("../../managerExample/DNA/ico.top","../../managerExample/DNA/ico.dat","DNA"); // Loading the top and the dat files
    cout <<ico.shiftbox({0,0,0})<<endl; // Bringing everything to 1st coordinate
    cout << ico.testBoxOverloaded()<<endl; // Checking overload from previous action.
    
    return 0;
}