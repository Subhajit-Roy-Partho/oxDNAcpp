#include "def.cpp"

int main(int argNum, char *argv[]) {

  Analysis crystal("output.top", "test4flat.dat", "crystal");
  crystal.correctA("./4crystals/newTest.dat");
  crystal.writeConfig("inputNew.dat");
  return 0;
}