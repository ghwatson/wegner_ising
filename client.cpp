/*
 * File: client.cpp
 * Author: ghwatson
 * Date: 25/05/2015
 */
#include <iostream>
#include "mc.hpp"
#include <boost/multi_array.hpp>
#include <kernel.hpp>
#include <fstream>

//General TODO list for project
//TODO: make parallelizable for use with job submission script
//TODO: pull together some graphing scripts (matplotlib?)
//rc file and matplotlib stuff from Chris!
//TODO: numerically verify Cv (20x20x20?)
//TODO: numerically verify Tc
//TODO: verify area+volume laws
//TODO: verify effects of disorder
//TODO: ultimately reproduce all the plots
//TODO: Binder cumulants
//TODO: single loop simulation 

using namespace std;

//Calculate the heat capacity for a variety of parameters.
int cv(){
  double Tf_start = 1.; //range of final temperatures
  double Tf_stop = 2.;
  double T_i = 10.; //start system here for each quench
  double step = 0.1;

  double e_i = 0;
  double e_f = 0.5;
  double step_e = 0.1;

  ofstream fout("cv.txt",ios::out);
  if (!fout){
    cout << "Failed to open stream." << endl;
    return 1;
  }

  //create the MC sim and the measurement kernel to insert into sim.
  WegnerMC sim = WegnerMC(e_i,T_i);
  KernelPipe* pipe = new KernelPipe();
  typedef void (KernelPipe::*kernel_ptr)(WegnerMC* sim); 
  kernel_ptr kernel = &KernelPipe::measure_Cv_data;

  //look at a range of disorders
  for (double e = e_i; e < e_f; e+=step_e){
    sim.set_e(e);
    cout << "for disorder " << e << endl;
    //we study a series of independent quenches (this first loop can be parallelized).
    //from T_i to Tf.
    for (double Tf = Tf_start; Tf < Tf_stop; Tf+=step){
      
      cout << "initializing to " << T_i << endl;
      sim.initialize(T_i); //initialize to random distribution at high temperature
      int eq = 10*10*10*10;
      sim.equilibrate(eq);

      //for each quench, we step from T_i to Tf slowly via 100 steps
      cout << "quenching to " << Tf << endl;
      double dT = (T_i - Tf) / 100.;
      for (double T = T_i; T >= Tf; T-=dT){ 
        sim.evolve(T,10*10); //evolve for 10^2
      }

      cout << "Done quench to " << Tf << ". Begin equilibration." << endl;

      //Now at Tf
      sim.equilibrate(eq); //equilibrate at desired temp

      cout << "Done eq. Begin measurement." << endl;

      int mc_steps = 10*10*10*10;
      sim.evolve(Tf, mc_steps, pipe, kernel); //collect data at desired temp

      cout << "Done measurement. Calc Cv" << endl;

      //Calculating Cv
      int N = 3*L*L*L;
      double Esq_avg = pipe->Esq_sum/mc_steps;
      double E_avg = pipe->E_sum/mc_steps;
      pipe->clean_data(); //re-initalize E, Esq
      double cv = (Esq_avg - E_avg*E_avg)/(Tf*Tf);
      cv = cv / N;
      fout << Tf << "      " << cv << "       " << e << endl;
    }//Tf
  }//e
  
  delete pipe;
  fout.close();

  return 0;
}

void wilson(){
  //wilson measurements
  
  WegnerMC sim = WegnerMC(10,0);

  //quench down to 1, measure -log(<W>).
  
}


void history(){
  //
}

int main(){
  //TODO: Put an argument if-then structure here eventually to choose scripts.
  int exit =  cv();
  std::cout << "done!" << std::endl;
  return exit;
}
