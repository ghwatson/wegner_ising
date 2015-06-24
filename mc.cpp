
#include <mc.hpp>
#include <math.h>
/*
 * bc -> unused OOO
 * sites -> unused EEE
 * edges -> used EEO
 * fc -> used EOO
 */
WegnerMC::WegnerMC(double e, double T){

  //Set disorder strength.
  m_e = e;
  m_T = T;

  //Allocate data structure.
  m_lattice.resize(boost::extents[2*L][2*L][2*L]);


  m_plaqs[0].resize( boost::extents[L][L][L]);
  m_plaqs[1].resize( boost::extents[L][L][L]);
  m_plaqs[2].resize( boost::extents[L][L][L]);

  range evens = range(0,2*L,2);
  range odds = range(1,2*L,2);

  //Edge Subviews (for 3D only)
  m_spins[0] = new view_3t(m_lattice[ boost::indices[odds][evens][evens] ]);
  m_spins[1] = new view_3t(m_lattice[ boost::indices[evens][odds][evens] ]);
  m_spins[2] = new view_3t(m_lattice[ boost::indices[evens][evens][odds] ]);

  //Store the disorder values on the plaquettes ( [normal][x][y][z] where x,y,z < N )
  m_disorders[0] = new view_3t(m_lattice[ boost::indices[evens][odds][odds] ]);
  m_disorders[1] = new view_3t(m_lattice[ boost::indices[odds][evens][odds] ]);
  m_disorders[2] = new view_3t(m_lattice[ boost::indices[odds][odds][evens] ]);

  //Initialize to ground state (all spins = 1)
  //TODO: figure out how to iterate subarrays (one-liner)
  for (int orient = 0; orient < n_dims; orient++)
    for (int i = 0; i < (int) m_spins[orient]->shape()[0];i++)
      for (int j = 0; j < (int) m_spins[orient]->shape()[1];j++)
        for (int k = 0; k < (int) m_spins[orient]->shape()[2];k++)
          (*m_spins[orient])[i][j][k] = 1;

  //Initialize plaquette disorders (randomly +/- 1, strength specified in e)
  for (int normal= 0; normal < n_dims; normal++)
    for (int i = 0; i < (int) m_disorders[normal]->shape()[0];i++)
      for (int j = 0; j < (int) m_disorders[normal]->shape()[1];j++)
        for (int k = 0; k < (int) m_disorders[normal]->shape()[2];k++)
          (*m_disorders[normal])[i][j][k] = ( (m_rgen.rand()>=0.5) ? +1 : -1 );

  //Initialize the plaquette products
  update_plaqs();

  //Get g.s. energy for future reference
  m_E0 = calc_E();

}//WegnerMC

WegnerMC::~WegnerMC(){
  //deallocate dynamically assigned memory.
  for (int i = 0; i < n_dims; i++){
    delete m_disorders[i];
    delete m_spins[i];
  }
}




// PUBLIC FUNCTIONS
   
  
  
  
  
void WegnerMC::evolve(double T, int mc_steps, void (*kernel)(WegnerMC*)){
  //do MC steps
  for (int mc = 0; mc < mc_steps; mc++){
    //do 3N^3 updates (=number of spins)
    for (int updates = 0; updates < 3*L*L*L; updates++){
      //pick a random spin.
      //random coordinates and orientation
      int orientation = m_rgen.randInt(2);
      int x = m_rgen.randInt(L);
      int y = m_rgen.randInt(L);
      int z = m_rgen.randInt(L);

      double expdelta = calc_dE(orientation, x, y, z, T);

      double r = m_rgen.rand();

      if (r > expdelta){
        int new_val = -(*m_spins[orientation])[x][y][z];
        (*m_spins[orientation])[x][y][z] = -new_val; //flip the spin
      }
    }
    //TODO: measurement kernel here
    //the kernel is a function which takes in a list of functions
    //the user should pass in pointers to the functions.
    kernel(this);
    
    //to measure Cv, calc E, E^2, and sum. after all loops done, div by N.
    

    
  }
  //update m_T
  m_T = T;
}//evolve

//void null_func(WegnerMC* sim);
void WegnerMC::equilibrate(int steps){
  //TODO: implement autocorrelation procedure to get this steps
  void null_func(WegnerMC* sim); //TODO: remove this prototype
  evolve(m_T, steps, null_func);
  std::cout << "afterwards" << std::endl;
}//equilibrate

void WegnerMC::initialize(double T_high){
  //Get random configuration
  for (int orient = 0; orient < n_dims; orient++)
    for (int i = 0; i < (int) m_disorders[orient]->shape()[0];i++)
      for (int j = 0; j < (int) m_disorders[orient]->shape()[1];j++)
        for (int k = 0; k < (int) m_disorders[orient]->shape()[2];k++)
          (*m_spins[orient])[i][j][k] = ( (m_rgen.rand()>=0.5) ? +1 : -1 );
  
  set_T(T_high);
  int steps = static_cast<int>(pow(10,4));
  //TODO: remove this
  steps = 20;
  equilibrate(steps);
}//initialize




// PRIVATE FUNS ----------





double WegnerMC::calc_E(){
  // Sum the terms in the Hamiltonian
  double E = 0;

  //update plaquettes
  update_plaqs();

  //sum all the plaquettes with their appropriate disorders.
  for (int orient = 0; orient < n_dims; orient++){
    for (int i = 0; i < L; i++){
      for (int j = 0; j < L;j++){
        for (int k = 0; k < L;k++){
          double ep = m_e*(*m_disorders[orient])[i][j][k];
          E -= (1 + ep)*m_plaqs[orient][i][j][k];
        }
      }
    }
  }
  return E;
}

double WegnerMC::calc_dE(int orientation, int x, int y, int z, int T){
  double expdelta;
  double E0 = 0;
  
  //loop over incident plaquettes by picking out vectors orthogonal to orientation
  for (int span_dir((orientation + 1) % 3); span_dir != orientation; span_dir = (span_dir + 1) % 3){

    //Plaquette spanned by orientation and span_dir
    int* spins = get_plaq_spins(orientation,span_dir,x,y,z);

    //Get the disorder term.
    int d_id[3] = {2*x, 2*y, 2*z};
    d_id[orientation] = (d_id[orientation] + 1) % L;
    d_id[span_dir] = (d_id[span_dir]+1) % L;
    boost::array<array_3t::index,3> id = {*d_id};
    double u = m_lattice(id);

    E0 += -(1+m_e*u)*spins[0]*spins[1]*spins[2]*spins[3];

    //Repeat for the negative plaquette
    int neg_site[3] = {x,y,z};
    neg_site[span_dir] = (L + neg_site[span_dir] - 1) % L;

    //This code is a repeat of the above, but using the negative indices instead
    //of x,y,z in order to access the appropriate plaquette.
    int* neg_spins = get_plaq_spins(orientation,span_dir,neg_site[0],neg_site[1],neg_site[2]);

    int neg_d_id[3] = {2*neg_site[0], 2*neg_site[1], 2*neg_site[2]};
    neg_d_id[orientation] = (neg_d_id[orientation] + 1) % L;
    neg_d_id[span_dir] = (neg_d_id[span_dir]+1) % L;
    boost::array<array_3t::index,3> neg_id = {*neg_d_id};
    double neg_u = m_lattice(neg_id);

    E0 += -(1+m_e*neg_u)*neg_spins[0]*neg_spins[1]*neg_spins[2]*neg_spins[3];

  }

  double E1 = -E0; //flipping the spin will flip the local energy.

  if (m_T >  0){
    expdelta = exp(-(E1-E0)/T);
  }
  else if (m_T == 0){
    if (E1 < E0){
      expdelta = 0;
    }
    else{
      expdelta = 1;
    }
  }
  else{
    std::cout << "Error with temperature value." << std::endl;
    return -1;
  }
  
  return expdelta;
}//calc_dE

double WegnerMC::calc_wilson(){
  //Calculate wilson loop.
  double wilson = 0;
  return wilson;
}

int* WegnerMC::get_plaq_spins(int orient, int span_dir, int x, int y, int z){
    //Positive plaquette spanned by orientation and span_dir.
    int* plaq_idx[4];
    plaq_idx[0] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[0][orient] = (plaq_idx[0][orient] + 1) % L;
    plaq_idx[1] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[1][span_dir] = (plaq_idx[1][span_dir]+1) % L;
    plaq_idx[1][orient] = (plaq_idx[1][orient] + 2) % L;
    plaq_idx[2] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[2][span_dir] = (plaq_idx[2][span_dir] + 2) % L;
    plaq_idx[3] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[3][span_dir] = (plaq_idx[3][span_dir] + 1) % L;


    //Get the spins on the plaquette
    int spins[4];
    for (int i = 0; i < 4; i++){
      boost::array<array_3t::index,3> id = {*plaq_idx[i]};
      spins[i] = m_lattice(id);
    }

    //TODO: replace plaq_idx with a multi_array.
    
    //Cleanup
    delete[] plaq_idx[0];
    delete[] plaq_idx[1];
    delete[] plaq_idx[2];
    delete[] plaq_idx[3];

    return spins;
}

void WegnerMC::update_plaqs(){
  //iterate across each site, and update the plaquettes
  for (int x = 0; x < L; x++){
  for (int y = 0; y < L; y++){
  for (int z = 0; z < L; z++){
  for (int orientation = 0; orientation < n_dims; orientation++){

    //For each orientation, pick a forward-permuted direction in the set {x,y,z}.
    //This covers all the positive plaquettes at a given site.
    int span_dir = (orientation + 1) % 3;

    int* spins = get_plaq_spins(orientation,span_dir,x,y,z);

    //Update
    int normal = (span_dir + 1) % 3;
    m_plaqs[normal][x][y][z] = spins[0]*spins[1]*spins[2]*spins[3];
   
  }//orientation
  }//z
  }//y
  }//x
}//update_plaqs


//Extras



//TODO:find a cleaner way to implement this null
  void null_func(WegnerMC* sim){
  }
