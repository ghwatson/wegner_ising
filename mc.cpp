
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

  range evens = range(0,2*L,2);
  range odds = range(1,2*L,2);

  //Edge Subviews (for 3D only)
  m_spins[0] = new view_3t(m_lattice[ boost::indices[odds][evens][evens] ]);
  m_spins[1] = new view_3t(m_lattice[ boost::indices[evens][odds][evens] ]);
  m_spins[2] = new view_3t(m_lattice[ boost::indices[evens][evens][odds] ]);
  m_test[0] = new view_3t(m_lattice[ boost::indices[odds][evens][evens] ]);

  //Store the disorder values on the plaquettes ( [normal][x][y][z] where x,y,z < N )
  m_disorders[0] = new view_3t(m_lattice[ boost::indices[evens][odds][odds] ]);
  m_disorders[1] = new view_3t(m_lattice[ boost::indices[odds][evens][odds] ]);
  m_disorders[2] = new view_3t(m_lattice[ boost::indices[odds][odds][evens] ]);

  m_plaqs[0].resize( boost::extents[2*L][2*L][2*L]);
  m_plaqs[1].resize( boost::extents[2*L][2*L][2*L]);
  m_plaqs[2].resize( boost::extents[2*L][2*L][2*L]);

  //Initialize to ground state (all spins = 1)
  std::fill(m_spins[0]->origin(), m_spins[0]->origin() + m_spins[0]->num_elements(), 3);
  std::fill(m_spins[1]->origin(), m_spins[1]->origin() + m_spins[1]->num_elements(), 3);
  std::fill(m_spins[2]->origin(), m_spins[2]->origin() + m_spins[2]->num_elements(), 3);
  //std::fill(m_test[0]->origin(), m_test[0]->origin() + m_test[0]->num_elements(), 333);
  std::cout << m_spins[0]->num_elements() << std::endl;
  std::cout << (*m_spins[1])[1][2][3] << std::endl;
  std::cout << (*m_test[0])[0][0][0] << std::endl;

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


  //TODO: DEBUG
  std::cout << "test value is: " << std::endl;
  (*m_spins[2])[0][0][0] = 23;
  std::cout << (*m_spins[2])[0][0][0] << std::endl;
  std::cout << (*m_spins[2])[0][1][2] << std::endl;
  std::cout << (m_spins[2]->origin()) << std::endl;
  std::cout << (m_lattice)[0][1][2] << std::endl;
  
}//WegnerMC

WegnerMC::~WegnerMC(){
  //deallocate dynamically assigned memory.
  for (int i = 0; i < n_dims; i++){
    //delete[] m_disorders[i];
    //delete[] m_spins[i];
  }
  std::cout << "balls" << std::endl;
}

// MEMBER FUNS ----------

double WegnerMC::calc_E(){
  // Sum the terms in the Hamiltonian
  double E = 0;

  
  return E;
}
//TODO: consider combinining this with evolve to enable the ability to
//perform local plaquette value updates after acceptance of a MC step.
//or pass the computed plaquette data out to the evolver.
double WegnerMC::calc_dE(int orientation, int x, int y, int z){
  double expdelta;
  double E0 = 0;
  
  //loop over incident plaquettes by picking out vectors orthogonal to orientation
  for (int span_dir((orientation + 1) % 3); span_dir != orientation; span_dir = (span_dir + 1) % 3){

    //Positive plaquette spanned by orientation and span_dir.
    int* plaq_idx[4];
    plaq_idx[0] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[0][orientation] = (plaq_idx[0][orientation] + 1) % L;
    plaq_idx[1] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[1][span_dir] = (plaq_idx[1][span_dir]+1) % L;
    plaq_idx[1][orientation] = (plaq_idx[1][orientation] + 2) % L;
    plaq_idx[2] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[2][span_dir] = (plaq_idx[2][span_dir] + 2) % L;
    plaq_idx[3] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[3][span_dir] = (plaq_idx[3][span_dir] + 1) % L;

    int d_id[3] = {2*x, 2*y, 2*z};
    d_id[orientation] = (d_id[orientation] + 1) % L;
    d_id[span_dir] = (d_id[span_dir]+1) % L;

    //Get the spins on the plaquette
    int spins[4];
    for (int i = 0; i < 4; i++){
      boost::array<array_3t::index,3> id = {*plaq_idx[i]};
      spins[i] = m_lattice(id);
    }

    //Get the disorder term.
    boost::array<array_3t::index,3> id = {*d_id};
    double u = m_lattice(id);

    E0 += -(1+m_e*u)*spins[0]*spins[1]*spins[2]*spins[3];
  }

  double E1 = -E0; //flipping the spin will flip the local energy.

  if (m_T >  0){
    expdelta = exp(-(E1-E0)/m_T);
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

void WegnerMC::evolve(double Tf, double dT){
  //Perform standard MC
  //get num_steps
  
  for (m_T = fabs(m_T - Tf); m_T < 0; m_T-=dT){
    //pick a random spin.
    //random coordinates and orientation
    int orientation = m_rgen.randInt(2);
    int x = m_rgen.randInt(L);
    int y = m_rgen.randInt(L);
    int z = m_rgen.randInt(L);

    double expdelta = calc_dE(orientation, x, y, z);

    double r = m_rgen.rand();

    if (r > expdelta){
      //int new_val = -(*m_spins[orientation])[x][y][z];
      //m_spins[orientation][x][y][z] = -new_val; //flip the spin
    }
    

    

    
    
    
    
    //get energy change for that spin
    //compare boltzmann to random number.
    //accept/refuse
    
    
    
    
    //perform a flip
    //accept/reject
    
  }
  //for each step, step in config space
  //lower energy?
  //if not then accept the change with boltzmann(Tf)
  //
  //Alter the data structure
  //Save all the data
  //User must choose observables from a list that will be output.
}

void WegnerMC::equilibrate(int steps){
  //TODO: implement autocorrelation procedure to get this steps
  evolve(m_T, steps);
}

////Convert 3D coordinates to 1D list index
//double WegnerMC::index_flatten(int x, int y, int z){
  //return 0;
//}

////Convert 1D list index to a 3D coordinate
//std::vector<int> WegnerMC::index_unflatten(int id){
  //std::vector<int> a;
  //return a;
//}

void WegnerMC::initialize(double T_high){
  //Get random configuration
  for (int orient = 0; orient < n_dims; orient++){
    std::fill(
    m_spins[orient]->origin(),
    m_spins[orient]->origin() + m_spins[orient]->num_elements(),
    (m_rgen.rand()>=0.5) ? +1 : -1);
  }

  set_T(T_high);
  int steps = static_cast<int>(pow(10,4));
  equilibrate(steps);
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

    //These are the indices for the spins on the plaquette
    
    int* plaq_idx[4];
    plaq_idx[0] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[0][orientation] = (plaq_idx[0][orientation] + 1) % L;
    plaq_idx[1] = new int[3] {2*x, 2*y, 2*z};
    plaq_idx[1][span_dir] = (plaq_idx[1][span_dir]+1) % L;
    plaq_idx[1][orientation] = (plaq_idx[1][orientation] + 2) % L;
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

    //Update
    int normal = (span_dir + 1) % 3;
    m_plaqs[normal][x][y][z] = spins[0]*spins[1]*spins[2]*spins[3];

    //Cleanup
    delete[] plaq_idx[0];
    delete[] plaq_idx[1];
    delete[] plaq_idx[2];
    delete[] plaq_idx[3];

    //delete[] plaq_idx;
  }//orientation
  }//z
  }//y
  }//x
}//update_plaqs
