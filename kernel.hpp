/*
 * File: kernel.hpp
 * Author: ghwatson
 * Date: 24/06/2015
 */
//This file details the implementations of the various kernels used in the simulation.
//
//They work by providing a class paired with a function.  The function is meant to be passed
//as a pointer to the evolve routine. One of its arguments will take the "this" pointer while
//inside of WegnerMC::evolve. This will allow the kernel to call various calculations of
//WegnerMC, and, since its a a member of a class external to WegnerMC, it will be capable of
//piping any sort of data outside of WegnerMC without needing change or overload the WegnerMC
//function.
#ifndef KERNEL_HPP
#define KERNEL_HPP

class WegnerMC; //forward declaration for mutual class defns

#include <boost/multi_array.hpp>

class KernelPipe{
  private:
    typedef boost::multi_array<double, 2> array_2t;
    typedef boost::multi_array<int, 1> array_1t;

  public:
    KernelPipe();
    void clean_data();

    //various kernels to use within Monte Carlo
    
    //--------CV kernel + variables
    void measure_Cv_data(WegnerMC* sim_pt);
    //TODO: add data to struct
    array_2t CvTe_data;
    double E_sum;
    double Esq_sum;


    //-------Wilson loop kernel
    void measure_wilson_data(WegnerMC* sim_pt);
    struct WilsonData{ //keep data organized in a struct
      //int wilson_sum;
      //int wilson_L; //TODO: could replace this with partial functions
      int Li;
      int Lf;
      array_1t loops;
    };
    WilsonData wilson_data;

    //-------history kernel
    
    //------hundred history kernel
    void measure_history_data(WegnerMC* sim_pt);
    struct HistoryData{
      array_2t history_data;
    };
};
#endif  /*KERNEL_HPP*/
