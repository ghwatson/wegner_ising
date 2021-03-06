/*
 * File: mc.hpp
 * Author: ghwatson
 * Date: 25/05/2015
 */

#ifndef WEGNER_MC_HPP
#define WEGNER_MC_HPP

#include <boost/multi_array.hpp>
#include "MersenneTwister.h"
#include <boost/foreach.hpp>

#include <typeinfo>

const int n_dims = 3;
//const int L = 8; //length in plaquettes/number of sites (PBC implies they're the same)

class KernelPipe; //forward declaration for mutual class defns

class WegnerMC{
  private:
    typedef boost::multi_array<int*, 6> array_6pt;
    //typedef boost::multi_array<int*, 5> array_5pt;
    typedef boost::multi_array<int, 4> array_4t;
    typedef boost::multi_array<int, 3> array_3t;
    typedef boost::multi_array<int, 2> array_2t;
    typedef boost::multi_array<int, 1> array_1t;
    typedef boost::multi_array_types::index_range range;
    //typedef array_5pt::array_view<4>::type view_4pt;
    typedef array_3t::array_view<3>::type view_3t;
    typedef array_2t::array_view<2>::type view_2t;
    
    MTRand m_rgen;
    int m_L;
    array_3t m_lattice;
    array_4t m_plaqs; //[normal][x][y][z]

    //Connectivity information
    array_6pt m_spin_nbs;    //[orient][x][y][z][plaq_id][spin_id]
                             //plaq_id is in {0..6} indexed by span_dir
                             //3..6 are the negative plaquettes

    //Provides an alternative way to access values via NxNxN site indices.
    view_3t* m_disorders[3]; //[normal][x][y][z] where x,y,z is incident site.
    view_3t* m_spins[3];  //[orientation][x][y][z]

    double m_E0; //g.s energy
    double m_T; //T
    double m_e; //disorder amount (will be +-m_e)

    double calc_dE(int orientation, int x, int y, int z, double T);
    void setup_connections();
    void update_plaqs();

  public:
    WegnerMC(int L, double e, double T);
    ~WegnerMC();
    double calc_E();
    void initialize(double T_high);
    void equilibrate(int steps);
    void evolve(double T, int mc_steps, KernelPipe* pipe, void (KernelPipe::*kernel)(WegnerMC*));
    void evolve(double T, int mc_steps);
    int get_plaq_val(int orient, int span_dir, int x, int y, int z);
    void set_T(double T);
    void set_e(double e);
};

#endif  /*WEGNER_MC_H*/
