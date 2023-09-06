#include "def.cpp"

int main(int argNum, char *argv[]) {

  Analysis crystal("output.top", "./4crystals/test.txt", "crystal");
  crystal.correctA("output.dat");
  crystal.writeConfig("./4crystals/newTest.dat");
  return 0;
}