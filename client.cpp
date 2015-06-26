/*
 * File: client.cpp
 * Author: ghwatson
 * Date: 25/05/2015
 */
#include <iostream>
#include "mc.hpp"
#include <boost/multi_array.hpp>
#include <kernel.hpp>

using namespace std;

//Calculate the heat capacity for a variety of parameters.
void cv(){
//TODO: is it possible to parallelize cv? or is sequential quenching needed
// to avoid metastable states?

  double Tf_start = 1; //range of final temperatures
  double Tf_stop = 2;
  double step = 0.05;
  double T_i = 10; //start system here for each quench

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
  
  KernelPipe* pipe = new KernelPipe();

  typedef void (KernelPipe::*kernel_ptr)(WegnerMC* sim); 
  kernel_ptr kernel = &KernelPipe::measure_Cv_data;

  //we study a series of independent quenches (this first loop can be
  //parallelized).
  for (double Tf = Tf_start; Tf < Tf_stop; Tf+=step){
    
    cout << "initializing to " << T_i << endl;
    sim.initialize(T_i); //initialize at high temperature

    //for each quench, we step from T_i to Tf slowly via 100 steps
    cout << "quenching to " << Tf << endl;
    double dT = (T_i - Tf) / 100.;
    for (double T = T_i; T > Tf; T-=dT){ 
      sim.evolve(Tf,10*10); //evolve for 10^2
      cout << "T is " << T << endl;
    }

    cout << "Done quench. Begin equilibration" << endl;

    //Now at Tf
    sim.equilibrate(10*10); //equilibrate at desired temp

    cout << "Done eq. Begin measurement." << endl;

    sim.evolve(Tf, 10*10, pipe, kernel); //collect data at desired temp

    cout << "Done measurement. Calc Cv" << endl;

    //Calculating Cv
    int N = 3*L*L*L;
    double Esq_avg = pipe->Esq_sum/100.;
    double E_avg = pipe->E_sum/100.;
    pipe->clean_data(); //re-initalize E, Esq
    double cv = (Esq_avg - E_avg*E_avg)/(Tf*Tf);
    cv = cv / N;
    cout << cv << " " << Tf << endl;
  }
  
  delete pipe;

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
  std::cout << "done!" << std::endl;
  return 0;
}
