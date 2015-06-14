#include <iostream>
#include "mc.hpp"
#include <boost/multi_array.hpp>

using namespace std;

//Calculate the heat capacity for a variety of parameters.
void cv(){
  WegnerMC sim = WegnerMC(0.1, 10);

  //start the temperature at 10

 
}

int main(){
  //TODO: Put an argument if-then structure here eventually to choose scripts.
  cv();
  return 0;
}
