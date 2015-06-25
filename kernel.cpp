/*
 * File: kernel.cpp
 * Author: ghwatson
 * Date: 25/06/2015
 */
#include <kernel.hpp>
#include <mc.hpp>

//typedef int (WegnerMC::*fnPtr)(int a, int b);

KernelPipe::KernelPipe(){
  //TODO: initialize here
  E_sum = 0;
  Esq_sum = 0;
}
void KernelPipe::measure_Cv_data(WegnerMC* sim_pt){
  double E_config = sim_pt->calc_E();
  E_sum += E_config;
  Esq_sum += E_config*E_config;
  std::cout << "E_sum = " << E_sum << std::endl;
  std::cout << "Esq_sum = " << Esq_sum << std::endl;
}
