/*
 * File: mc.hpp
 * Author: ghwatson
 * Date: 25/05/2015
 */

#ifndef WEGNER_MC_H_
#define WEGNER_MC_H_

#include <boost/multi_array.hpp>
#include "MersenneTwister.h"
#include <boost/foreach.hpp>

#include <typeinfo>

const int n_dims = 3;
const int L = 10; //length in plaquettes/number of sites (PBC implies they're the same)

class WegnerMC{
  private:
    typedef boost::multi_array<int, 3> array_3t;
    typedef boost::multi_array<int, 1> array_1t;
    typedef boost::multi_array_types::index_range range;
    typedef array_3t::array_view<3>::type view_3t;
    
    MTRand m_rgen = MTRand();
    array_3t m_lattice;
    array_3t m_plaqs[3];

    
    //TODO: replace with 4d arrays to avoid pointy stuff, and * in calling values.
    //Provides an alternative way to access values via NxNxN site indices.
    view_3t* m_disorders[3]; //[normal][x][y][z] where x,y,z is incident site.
    view_3t* m_spins[3];  //[orientation][x][y][z]

    double m_T; //T
    double m_E0; //g.s energy
    double m_e; //disorder amount (will be +-m_e)

    //void calc_plaq(int normal, int x, int y, int z);
    double calc_dE(int orientation, int x, int y, int z, int T);
    double calc_wilson();

    int* get_plaq_spins(int orient, int span_dir, int x, int y, int z);

    void update_plaqs();

  public:
    WegnerMC(double e, double T);
    ~WegnerMC();
    double calc_E();
    void initialize(double T_high);
    void equilibrate(int steps);
    //void evolve(double Tf, int steps, void (*measure)());
    void evolve(double T, int mc_steps);
    void evolve(double T, int mc_steps, void (*kernel)(WegnerMC*));
    void set_T(double T){m_T = T;};
    void set_e(double e){m_e = e;};
};

#endif  /*WEGNER_MC_H_*/
