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
const int N = 20; //TODO: template these in.

//boost boilerplate
// Using boost arrays to make it easier to generalize to higher dimensions.
//typedef boost::multi_array<int, n_dims> array_nt;
//typedef boost::multi_array<int, 3> array_3t;
//typedef boost::multi_array<int, 2> array_2t;
//typedef boost::multi_array<int, 1> array_1t;
//typedef boost::multi_array_types::index_range range_t;
//typedef boost::multi_array<int, 1>::const_iterator c_iter;
//typedef array_nt::array_view<n_dims> subView;

// Tools for multi-arrays.
//array_1t::extent_gen extents;
//array_2t::index_gen indices;

//TODO: To generalize this to hypertoric codes, will need to swap out data
//structure for a combinatorial map or linear complex data structure (see CGAL)
//template <int n_dims>
class WegnerMC{
  private:
    //typedef boost::multi_array<int, n_dims> array_nt;
    typedef boost::multi_array<int, 3> array_3t;
    typedef boost::multi_array_types::index_range range;
    typedef array_3t::array_view<3>::type view_3t;
    
    MTRand m_rgen = MTRand();
    array_3t m_lattice;

    //TODO: In future, a generic scheme would have some tuple paired with each lattice point
    //so that the face centers can store both disorder and product value.
    array_3t m_plaqs[3]; //plaquette products

    // Views for the plaquettes disorder terms and the spins
    //TODO: use dynamic arrays to move this junk to constructor?
    //view_3t m_disorders[3] = {
        //m_lattice(boost::indices[evens][odds][odds]),
        //m_lattice(boost::indices[odds][evens][odds]),
        //m_lattice(boost::indices[odds][odds][evens])};
    //view_3t m_spins[3] = {
        //m_lattice(boost::indices[odds][evens][evens]),
        //m_lattice(boost::indices[evens][odds][evens]),
        //m_lattice(boost::indices[evens][evens][odds])};
    
    //view_3t* testview;
    view_3t** m_disorders = new view_3t*[3];
    view_3t** m_spins = new view_3t*[3];

    //Neighbourhood + Incidence relationships
    array_1t m_spin_to_plaqs;
    array_1t m_plaq_to_spins;

    double m_Th; //= 100;
    double m_E0; //g.s energy

    //TODO: Change structure to utilize local updates to get dE.
    double calc_E();
    double calc_plaq(int normal, int x, int y, int z);
    double calc_dE();
    double calc_Eflucs();
    double calc_Cv();
    double calc_wilson();

    //Index conversion
    double index_flatten(int x, int y, int z);
    std::vector<int, 3> index_unflatten(int id);

    void update_plaqs();
    void update_plaqs(int x, int y, int z);
    void update_plaqs(int orientation, int x, int y, int z);
 
  public:
    WegnerMC(double e);
    void initialize(double T_high);
    void evolve(double T, int steps);
    void set_T_high(double Th){m_Th = Th;};

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
//TODO: How would you write generic multi_array code without specifying dimension
//template <int n_dims>
//WegnerMC<n_dims>::WegnerMC(int N, double e){
WegnerMC::WegnerMC(double e){

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
  std::fill(
      m_disorders[0]->origin(),
      m_disorders[0]->origin() + m_disorders[0]->num_elements(),
      (m_rgen.rand()>=0.5) ? +1 : -1);
  std::fill(
      m_disorders[1]->origin(),
      m_disorders[1]->origin() + m_disorders[1]->num_elements(),
      (m_rgen.rand()>=0.5) ? +1 : -1);
  std::fill(
      m_disorders[2]->origin(),
      m_disorders[2]->origin() + m_disorders[2]->num_elements(),
      (m_rgen.rand()>=0.5) ? +1 : -1);

  //Initialize the plaquette products
  update_plaqs();


  //Get g.s. energy for future reference
  m_E0 = calc_E();
  
}//WegnerMC(N,e)

//TODO: make a function that updates plaquette values, given a local update.


//template <int n_dims>
//double WegnerMC<n_dims>::calc_E(){
double WegnerMC::calc_E(){
  // Sum the terms in the Hamiltonian
  double E = 0;

  
  return E;
}
//template <int n_dims>
//double WegnerMC<n_dims>::calc_dE(){
double WegnerMC::calc_dE(){
  //Get energy given previous energy step.
  //pass in an array to allocate to?
  double dE = 0;
  return dE;
}
//template <int n_dims>
//double WegnerMC<n_dims>::calc_Eflucs(){
double WegnerMC::calc_Eflucs(){
  //stuff
  //<E^2> - <E>^2
  double Eflucs = 0;
  return Eflucs;
}
//template <int n_dims>
//double WegnerMC<n_dims>::calc_Cv(){
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


//Convert 3D coordinates to 1D list index
double WegnerMC::index_flatten(int x, int y, int z){
  return 0;
}

//Convert 1D list index to a 3D coordinate
std::vector<int, 3> index_unflatten(int id){
  return 0;
}


//TODO: Is this function useful?
void WegnerMC::update_plaqs(){
}
void WegnerMC::update_plaqs(int x, int y, int z){
  //For accessing via 2Nx2Nx2N lattice
  int product = 1;
  int n_nbs = 2*n_dims;
  int plaq1 = m_lattice[x][y][z]
      *m_lattice[x][y][z+2]
      *m_lattice[x][y+1][z+1]
      *m_lattice[x][y-1][z+1];

  

  //int plaq2 =
  //for (int i = 0; i < n_nbs; i++){
    //product *= m_lattice[
  //}


}
void WegnerMC::update_plaqs(int orientation, int x, int y, int z){
  //Given the spin site, update all neighbouring plaquettes.
  //(*m_plaqs[0])(idx) = 
  //int product = 1;
  //int n_nbs = 2*n_dims;
  
  //iterate over nbs_plaqs
  
  //STEPS:
  //1. access an fc site by picking a direction perpendicular to the variable direction
  //2. from the fc site, step out in all directions in the plane spanned by "direction"
  //    and fc direction.
  //3. Repeat for other directions perpendicular to "direction"

  //for (int i = 0; i < n_nbs; i++){
    //product *= m_spins[direction]
  
  //Using the direction, figure out the relative coordinates of the fc faces.
  for (int dir = 0; dir < n_dims; dir++){

    if (dir == orientation){
      continue;
    }


  }

    
    //positive direction
    int fc1[3] = {x+1, y, z};
    int fc2[3] = {x, y, z};
   
    //negative direction
  
  //Convert to lattice index.
  //Then just shift from these fcs to get the neighbouring spins!
  
}

//template <int n_dims>
//void WegnerMC<n_dims>::Initialize(double T_high){
void WegnerMC::initialize(double T_high){
  //Get random configuration

  //Evolve until equilibrated
}

//template <int n_dims>
//void WegnerMC<n_dims>::Evolve(double T, int steps){
void WegnerMC::evolve(double T, int steps){
  //Perform standard MC
  //Alter the data structure
  //Save all the data
  //User must choose observables from a list that will be output.
}

//void convert_index(){
  //// change between 1d and 3d index using mod.
  //// is this necessary?
//}


//int main() {

	//cout << "!!!Hello World!!!" << endl; 
  ////WegnerMC<3> sim(10,0.1);

	//return 0;
//}


//typedef boost::multi_array<int, 3> array_3t;

//TODO: wrap into class later?
//TODO: later generalize to templated array without dimension defined.
//Subview functions for accessing the data structure

//typedef array_nt::array_view<3> subView;
//
//class WegnerData{
//private:
//public:
//  array_nt m_lattice;
//  subViiew m_edges;
//  subViiew m_sites;
//  subViiew m_fc;
//  subViiew m_bc;
//
//  WegnerData(int N){
//    //Allocate space.
//    m_lattice.resize(boost::extents[2*N][2*N][2*N]);
//
//    //Create the subviews.
////
//// array_view dims:
//// [base,stride,bound)
//// [0,1,2), [1,1,3), [0,2,4)
////
//
////  typedef array_3t::index_range range;
////  m_edges = array_3t[boost::indices[range(0,2)][range(1,3)][range(0,4,2)]];
////
////  for (array::index i = 0; i != 2; ++i)
////    for (array::index j = 0; j != 2; ++j)
////      for (array::index k = 0; k != 2; ++k)
////        assert(myview[i][j][k] == myarray[i][j+1][k*2]);
//
//  }

//};
#endif  /*WEGNER_MC_H_*/
