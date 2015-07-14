/*
 * File: kernel.cpp
 * Author: ghwatson
 * Date: 25/06/2015
 */
#include <kernel.hpp>
#include <mc.hpp>

using namespace std;

KernelPipe::KernelPipe(){
  clean_data();
}

//re-initializes all data
void KernelPipe::clean_data(){
  E_sum = 0;
  Esq_sum = 0;
  //cv_data = CvData();

  wilson_data = WilsonData();
}
void KernelPipe::measure_Cv_data(WegnerMC* sim_pt){
  double E_config = sim_pt->calc_E();
  E_sum += E_config;
  Esq_sum += E_config*E_config;
  //cout << "econfig = " << E_config << endl;
  //std::cout << "E_sum = " << E_sum << std::endl;
  //std::cout << "Esq_sum = " << Esq_sum << std::endl;
}
void KernelPipe::measure_wilson_data(WegnerMC* sim_pt){
  //measure square loop of size L in the x-y plane at z = 0
  int wilson_loop = 1;
  for (int L = wilson_data.Li; L <= wilson_data.Lf; L++){
    for (int x = 0; x < L; x++){
      for (int y = 0; y < L; y++){
        //get the plaquette in the x,y plane
        //cout << wilson_loop << endl;
        wilson_loop *= sim_pt->get_plaq_val(0,1,x,y,0);
        //cout << " --------" << endl;
        //cout <<  wilson_data.loops[L] << endl;
        wilson_data.loops[L-1] += wilson_loop;
        //cout <<  wilson_data.loops[L-1] << endl;
        if (wilson_loop == 0){
          cout << "TROUBLE" << endl;
          getchar();
        }
      }
    }
  }
}

void KernelPipe::measure_history_data(WegnerMC* sim_pt){
  //get the dE values, but dont sum, just pop onto a vector.

}
