#include <iostream>
#include "mc.hpp"
#include <boost/multi_array.hpp>

using namespace std;
typedef boost::multi_array<int, 3> array_3t;

void myfunc(array_3t anarray){
  anarray.resize(boost::extents[3][2][1]);
}

class testclass{
  array_3t thisarray;
  public:
    testclass(int x, int y, int z);
};

testclass::testclass(int x,int y,int z){
  thisarray.resize(boost::extents[x][y][z]);
}

int main(){
  cout << "hello world" << endl;
  array_3t m(boost::extents[1][2][3]);

  myfunc(m);

  testclass myobj(1,3,2);
  
  WegnerMC sim = WegnerMC(0.2);

  cout << "success!" << endl;
  
  return 0;
}

