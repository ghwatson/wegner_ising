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
const int N = 20;

class WegnerMC{
  private:
    typedef boost::multi_array<int, 3> array_3t;
    typedef boost::multi_array<int, 1> array_1t;
    typedef boost::multi_array_types::index_range range;
    typedef array_3t::array_view<3>::type view_3t;
    
    MTRand m_rgen = MTRand();
    array_3t m_lattice;
    array_3t m_plaqs[3]; //plaquette products
    
    view_3t** m_disorders = new view_3t*[3];
    view_3t** m_spins = new view_3t*[3];

    //Neighbourhood + Incidence relationships
    //TODO: Implement these? (move logic from update_plaqs)
    //array_1t m_spin_to_plaqs;
    //array_1t m_plaq_to_spins;

    double m_Ti; //T_initial
    double m_E0; //g.s energy
    double m_e;

    //TODO: Change structure to utilize local updates to get dE.
    double calc_E();
    double calc_plaq(int normal, int x, int y, int z);
    double calc_dE();
    double calc_Eflucs();
    double calc_Cv();
    double calc_wilson();

    //Index conversion
    double index_flatten(int x, int y, int z);
    std::vector<int> index_unflatten(int id);

    void update_plaqs();
 
  public:
    WegnerMC(double e);
    ~WegnerMC();
    void initialize(double T_high);
    void evolve(double Tf, int steps);
    void set_T_high(double Ti){m_Ti = Ti;};
    void set_e(double e){m_e = e;};

    bool flag_plaq = 1; //Spin flips will locally update plaquette products
};




/*
 * IMPLEMENTATION
 */


/*
 * bc -> unused OOO
 * sites -> unused EEE
 * edges -> used EEO
 * fc -> used EOO
 */
WegnerMC::WegnerMC(double e){

  //Set disorder strength.
  m_e = e;

  //Allocate data structure.
  m_lattice.resize(boost::extents[2*N][2*N][2*N]);

  range evens = range(0,2*N,2);
  range odds = range(1,2*N,2);

  //Edge Subviews (for 3D only)
  m_spins[0] = new view_3t(m_lattice[ boost::indices[odds][evens][evens] ]);
  m_spins[1] = new view_3t(m_lattice[ boost::indices[evens][odds][evens] ]);
  m_spins[2] = new view_3t(m_lattice[ boost::indices[evens][evens][odds] ]);

  //Plaquettes Subviews (fc is characterized by OOE,OEO,OOE... m_plaqs[normal direction][x][y][z])
  m_disorders[0] = new view_3t(m_lattice[ boost::indices[evens][odds][odds] ]);
  m_disorders[1] = new view_3t(m_lattice[ boost::indices[odds][evens][odds] ]);
  m_disorders[2] = new view_3t(m_lattice[ boost::indices[odds][odds][evens] ]);

  //Initialize to ground state (all spins = 1)
  std::fill(m_spins[0]->origin(), m_spins[0]->origin() + m_spins[0]->num_elements(), 1);
  std::fill(m_spins[1]->origin(), m_spins[1]->origin() + m_spins[1]->num_elements(), 1);
  std::fill(m_spins[2]->origin(), m_spins[2]->origin() + m_spins[2]->num_elements(), 1);

  //Initialize plaquette disorders (randomly +/- 1, strength specified in e)
  for (int normal = 0 ; normal < n_dims; normal++){
    std::fill(
        m_disorders[normal]->origin(),
        m_disorders[normal]->origin() + m_disorders[normal]->num_elements(),
        (m_rgen.rand()>=0.5) ? +1 : -1);
  }
  
  //Initialize the plaquette products
  update_plaqs();

  //Get g.s. energy for future reference
  m_E0 = calc_E();
  
}//WegnerMC

WegnerMC::~WegnerMC(){
  //deallocate dynamically assigned memory.
  for (int i = 0; i < n_dims; i++){
    delete[] m_disorders[i];
    delete[] m_spins[i];
  }
}

double WegnerMC::calc_E(){
  // Sum the terms in the Hamiltonian
  double E = 0;

  
  return E;
}
double WegnerMC::calc_dE(){
  //Get energy given previous energy step.
  //pass in an array to allocate to?
  double dE = 0;
  return dE;
}
double WegnerMC::calc_Eflucs(){
  //stuff
  //<E^2> - <E>^2
  double Eflucs = 0;
  return Eflucs;
}
double WegnerMC::calc_Cv(){
  //stuff
  double Cv = 0;
  return Cv;
}
double WegnerMC::calc_wilson(){
  //Calculate wilson loop.
  double wilson = 0;
  return wilson;
}

void WegnerMC::evolve(double Tf, int steps){
  //Perform standard MC
  //Alter the data structure
  //Save all the data
  //User must choose observables from a list that will be output.
}

//Convert 3D coordinates to 1D list index
double WegnerMC::index_flatten(int x, int y, int z){
  return 0;
}

//Convert 1D list index to a 3D coordinate
std::vector<int> WegnerMC::index_unflatten(int id){
  std::vector<int> a;
  return a;
}

void WegnerMC::initialize(double T_high){
  //Get random configuration

  //Evolve until equilibrated
}

void WegnerMC::update_plaqs(){
  //iterate across each site, and update the plaquettes
  for (int x = 0; x < N; x++){
  for (int y = 0; y < N; y++){
  for (int z = 0; z < N; z++){
  for (int orientation = 0; orientation < n_dims; orientation++){
    //Work with the 2Nx2Nx2N index.
    int indices_2d[3] = {2*x,2*y,2*z};

    //For each orientation, pick a forward-permuted direction in the set {x,y,z}.
    //This covers all the positive plaquettes at a given site.
    int span_dir = (orientation + 1) % 3;

    //These are the indices for the spins on the plaquette
    int* plaq_idx[4];
    plaq_idx[0] = indices_2d;
    plaq_idx[0][orientation] = (plaq_idx[0][orientation] + 1) % N;
    plaq_idx[1] = indices_2d;
    plaq_idx[1][span_dir] = (plaq_idx[1][span_dir]+1) % N;
    plaq_idx[1][orientation] = (plaq_idx[1][orientation] + 2) % N;
    plaq_idx[2] = indices_2d;
    plaq_idx[2][span_dir] = (plaq_idx[2][span_dir] + 2) % N;
    plaq_idx[3] = indices_2d;
    plaq_idx[3][span_dir] = (plaq_idx[3][span_dir] + 1) % N;

    //Get the spins on the plaquette
    int spins[4];
    for (int i = 0; i < 4; i++){
      boost::array<array_3t::index,3> id = {*plaq_idx[i]};
      spins[i] = m_lattice(id);
    }

    //Update
    int normal = (span_dir + 1) % 3;
    m_plaqs[normal][x][y][z] = spins[0]*spins[1]*spins[2]*spins[3];

    if (x == 20 and y == 0 and z == 0 and orientation == 2){
      std::cout << "hey0" << std::endl;
    }
  }//orientation
  }//z
  }//y
  }//x
}//update_plaqs

#endif  /*WEGNER_MC_H_*/
