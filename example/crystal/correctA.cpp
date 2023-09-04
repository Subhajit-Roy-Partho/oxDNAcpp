#include "def.cpp"

int main(int argNum, char *argv[]) {

  Analysis crystal("output.top", "test.dat", "crystal");
  crystal.correctA("output.dat");
  crystal.writeConfig("newTest");
  return 0;
}