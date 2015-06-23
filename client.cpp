#include <iostream>
#include "mc.hpp"
#include <boost/multi_array.hpp>

using namespace std;

//Calculate the heat capacity for a variety of parameters.
void cv(){
//TODO: is it possible to parallelize cv? or is sequential quenching needed
// to avoid metastable states?

  double Tf_start = 1;
  double Tf_stop = 2;
  double step = 0.05;

  //start the temperature at 10, disorder 0
  WegnerMC sim = WegnerMC(0, 10);

  //for a range of Tfinals
  //start at 10
  //equilibrate
  //for a step of dT towards Tfinal
  //short eq of 10^2
  //10^4 eq at Tfinal
  //do measure over 10^4 steps
  //record the (T,Cv,e) value to file
  
  
  for (double Tf = Tf_start; Tf < Tf_stop; Tf+=step){
    sim.initialize(10);
  }

}

void wilson(){
  //wilson measurements
}

void history(){
  //
}

int main(){
  //TODO: Put an argument if-then structure here eventually to choose scripts.
  cv();
  return 0;
}
