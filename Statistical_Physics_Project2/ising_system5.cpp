#include "ising.h"

int main(int argc, char *argv[])
{
  
  if(argc != 3)
  {
    cout << "please supply the arguments L , n_steps" << endl;
    exit(0);
  }  

  // Number of particle aling a row
  int L = (int) strtol(argv[1], NULL, 10); 
  // number of thermalization to be done
  int n_steps = strtol(argv[2], NULL, 10);

  // Number of neigbohrs for each particle
  int n_br = 4;
  // Neigbohrs of each particle
  int ** nbr;
  // Spins of each particle
  int * spins;
  // allocating memory
  allocate_system(L, n_br, nbr, spins);
  // name of file to save data
  string name;

  

  /* Init and running modelation for a given seed for thermalization and
      different for spins */
  int i;
  double beta;
/*  cout << "----------------------------------------------------------" << "\n";
  cout << "Results for same thermalization different initials beta 0" << "\n";
  beta = 0;
  for(i = 0; i < 2; i ++){
    srand(i);  // choosing seed
    init_periodic_system(L, n_br, nbr, spins);

    name = "Periodic_data/alleatory1beta1seed1";
    moments(L, n_br, n_steps, beta, nbr, spins, name, 1);
  }

  // For up and down spins
  for(i = 0; i < (L * L ); i++)
    spins[i] = 1;
  
  moments(L, n_br, n_steps, beta, nbr, spins, name, 1);

  // For up and down spins
  for(i = 0; i < (L * L ); i++)
    spins[i] = -1;
  
  moments(L, n_br, n_steps, beta, nbr, spins, name, 1);
*/
  
  cout << "----------------------------------------------------------" << "\n";
  cout << "Results for same thermalization different initials beta 1" << "\n";

  beta = 1;
  init_periodic_system(L, n_br, nbr, spins);
  
  for(int j = 200; j < 201; ++j){ 
      for(i = 0; i < 2; i ++){
        srand(j);  // choosing seed
        init_periodic_system(L, n_br, nbr, spins);

        name = "Periodic_data/alleatory1beta1seed1";
        moments(L, n_br, n_steps, beta, nbr, spins, name, 1);
      }

  }

  /*// For up and down spins
  for(i = 0; i < (L * L ); i++)
    spins[i] = 1;
  
  moments(L, n_br, n_steps, beta, nbr, spins, name, 1);

  // For up and down spins
  for(i = 0; i < (L * L ); i++)
    spins[i] = -1;
  
  moments(L, n_br, n_steps, beta, nbr, spins, name, 1);*/


  /* Init and running modelation for a different seed for thermalization and
      different for spins */
/*
  beta = 0;
  cout << "----------------------------------------------------------" << "\n";
  cout << "Results for same thermalization process different initials" << "\n";
  for(i = 0; i < 2; i ++){
    srand(i);  // choosing seed
    init_periodic_system(L, n_br, nbr, spins);

    name = "Periodic_data/alleatory1beta1seed1";
    moments(L, n_br, n_steps, beta, nbr, spins, name, i);
  }

  // For up and down spins
  for(i = 0; i < (L * L ); i++)
    spins[i] = 1;
  
  moments(L, n_br, n_steps, beta, nbr, spins, name, 1);

  // For up and down spins
  for(i = 0; i < (L * L ); i++)
    spins[i] = -1;
  
  moments(L, n_br, n_steps, beta, nbr, spins, name, 1);*/



}