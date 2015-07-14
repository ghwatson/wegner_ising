/*
 * File: client.cpp
 * Author: ghwatson
 * Date: 25/05/2015
 */
#include <iostream>
#include <sstream>
#include <iomanip>
#include "mc.hpp"
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>
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
    return EXIT_FAILURE;
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
      double Esq_avg = pipe->Esq_sum/mc_steps;
      double E_avg = pipe->E_sum/mc_steps;
      pipe->clean_data(); //re-initalize E, Esq
      double cv = (Esq_avg - E_avg*E_avg)/(Tf*Tf);
      int N = 3*L*L*L;
      cv = cv / N;
      fout << Tf << "      " << cv << "       " << e << endl;
    }//Tf
  }//e
  
  delete pipe;
  fout.close();

  return EXIT_SUCCESS;
}

int job_cv(double e,double Tf){
  double T_i = 10.; //start system here for the  quench

  //Create file to output to for this quench
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(5);
  ss << "outputdata/e_" << e << "_T_" << Tf << ".txt";
  std::string file_name = ss.str();
  ofstream fout(file_name, ios::out);
  if (!fout){
    cout << "Failed to open stream." << endl;
    return EXIT_FAILURE;
  }

  //create the MC sim and the measurement kernel to insert into sim.
  WegnerMC sim = WegnerMC(e,T_i);
  KernelPipe* pipe = new KernelPipe();
  typedef void (KernelPipe::*kernel_ptr)(WegnerMC* sim); 
  kernel_ptr kernel = &KernelPipe::measure_Cv_data;

  sim.initialize(T_i); //initialize to random distribution at high temperature
  int eq = 10*10*10*10;
  sim.equilibrate(eq);

  //for each quench, we step from T_i to Tf slowly via 100 steps
  double dT = (T_i - Tf) / 100.;
  for (double T = T_i; T >= Tf; T-=dT){ 
    sim.evolve(T,10*10); //evolve for 10^2
  }

  //Now at Tf
  sim.equilibrate(eq); //equilibrate at desired temp

  int mc_steps = 10*10*10*10;
  sim.evolve(Tf, mc_steps, pipe, kernel); //collect data at desired temp

  //Calculating Cv
  double Esq_avg = pipe->Esq_sum/mc_steps;
  double E_avg = pipe->E_sum/mc_steps;
  pipe->clean_data(); //re-initalize E, Esq
  double cv = (Esq_avg - E_avg*E_avg)/(Tf*Tf);
  int N = 3*L*L*L;
  cv = cv / N;
  fout << Tf << "      " << cv << "       " << e << endl;
  
  delete pipe;
  fout.close();

  return EXIT_SUCCESS;
}

int wilson(){
  //TODO: finish and debug this.
  //wilson measurements for a clean system
  
  WegnerMC sim = WegnerMC(10,0);
  KernelPipe* pipe = new KernelPipe();
  typedef void (KernelPipe::*kernel_ptr)(WegnerMC* sim); 
  kernel_ptr kernel = &KernelPipe::measure_wilson_data;
  
  double eq = 10*10*10*10;

  ofstream fout("wilson.txt",ios::out);
  if (!fout){
    cout << "Failed to open stream." << endl;
    return EXIT_FAILURE;
  }

  double T_i = 10;
  sim.initialize(T_i); //set random state
  sim.equilibrate(eq); //equilibrate after random state
  //quench down to 1.5
  double T1 = 1.5;
  double dT = (T_i - T1) / 100.;
  for (double T = T_i; T >= T1; T-=dT){ 
    sim.evolve(T,10*10); //evolve for 10^2
  }
  sim.equilibrate(eq); //equilibrate at desired temp

  //measure wilson loops
  int Li = 1; int Lf = 5;
  pipe->wilson_data.Li = Li;
  pipe->wilson_data.Lf = Lf;
  pipe->wilson_data.loops.resize(boost::extents[Lf-Li + 1]);
  int mc_steps = 10*10*10*10;
  sim.evolve(T1, mc_steps, pipe, kernel); //collect data at desired temp
  //Now get the averages and write to file
  for (int L = Li; L <= Lf; L++){ //Note this loop only works for dL = 1
    int wilsonL_avg = pipe->wilson_data.loops[L-1] / double(mc_steps);
    fout << T1 << " " << L << " " << wilsonL_avg << endl;
  }
  pipe->clean_data();
  


  //now repeat wilson measurements, but at T = 1.
  double T2 = 1;
  dT = (T1 - T2) / 100.; //TODO: worried about the steps here?
  for (double T = T1; T >= T2; T-=dT){ 
    sim.evolve(T,10*10); //evolve for 10^2
  }
  sim.equilibrate(eq); //equilibrate at desired temp
  //measure wilson loops
  pipe->wilson_data.Li = Li;
  pipe->wilson_data.Lf = Lf;
  pipe->wilson_data.loops.resize(boost::extents[Lf-Li + 1]);
  sim.evolve(T2, mc_steps, pipe, kernel); //collect data at desired temp
  //Now get the averages and write to file
  for (int L = Li; L <= Lf; L++){ //Note this loop only works for dL = 1
    double wilsonL_avg = pipe->wilson_data.loops[L-1] / double(mc_steps);
    fout << T2 << " " << L << " " << wilsonL_avg << endl;
  }

  //quench down to 1.5, measure -log(<W>) for various sizes.
  //quench down to 1, measure -log(<W>).
  
  delete pipe;
  
  return EXIT_SUCCESS;
}

void individual_history(){
  //for a set of initial configurations 1-10. for now, choose random initial configurations.
  //for a range of e.
  //for a range of Tf to quench to.
  //measure dE at each MCstep. the step acts as time t.
  //record Tf, t, dE, e(k?), trial#
}

void hundred_history(){
  //do a 100 random configurations.
  for (int i = 0; i < 100; i++){
  }

  //for a range of e. note that the randomization is chosen differently each time!
  //TODO: is the randomization chosen differently elsewhere???
  //for a range of Tf to quench to.
  //measure dE at each MCstep. the step acts as time t.
  //record Tf, t, dE, e(k?), trial#
}

int main(int argc, const char* argv[]){

  // CV stuff
  double e = atof(argv[1]);
  double Tf = atof(argv[2]);

  int exit = job_cv(e,Tf);

  return exit;
}
