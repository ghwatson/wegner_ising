#include <mc.hpp>


//This file details the implementations of the various kernels used in the simulation.
//
//They work by providing a class paired with a function.  The function is meant to be passed as a pointer to the evolve routine. One of its arguments will take the "this" pointer while inside of WegnerMC::evolve. This will allow the kernel to call various calculations of WegnerMC, and, since its a a member of a class external to WegnerMC, it will be capable of piping any sort of data outside of WegnerMC without needing change or overload the WegnerMC function.



class KernelPipe{
  private:
    typedef boost::multi_array<double, 2> array_2t;


  public:
    KernelPipe();

    //--------CV kernel + variables
    //various kernels to use within Monte Carlo
    void measure_Cv_data(WegnerMC* sim_pt);
    void testkernel();
    array_2t CvTe_data;
    double E_sum = 0;
    double Esq_sum = 0;

};

KernelPipe::KernelPipe(){

  std::cout << "hello" << std::endl;
}

void KernelPipe::testkernel(){}

//typedef int (WegnerMC::*fnPtr)(int a, int b);
void KernelPipe::measure_Cv_data(WegnerMC::* sim_pt){
  std::cout << "hello" << std::endl;
  double E_config = sim_pt->calc_E();
  E_sum += E_config;
  Esq_sum += E_config*E_config;
  std::cout << E_sum << std::endl;
  std::cout << Esq_sum << std::endl;
  std::cout << "hello" << std::endl;
}
