#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>
#include "matrix.h"
#include "elasticity.h"
#include "math.h"
#include "lu.h"
#include "design.h"
#include "eigen.h"
#include "optimization.h"
#include "animation.h"

int main (int argc, char *argv[]) {

  if (argc < 2){
    printf("Usage: \n"
			"./project <results_filename> <animation_filename><\n" 
			"---------------------------- \n\n"
			"- result_filename is the output file to write the dimensions of the tuningFork \n "
			"- optional : animation_filename is the output file to save the animation (.mp4)\n "
      "\n");
		return -1;
  } 

  char results_filename[100];
  sprintf(results_filename, "results/%s",argv[1]);
  char * animation_filename;

  int animation = 0;
  if (argc == 3) {
    animation = 1;
    animation_filename = argv[2];
  }

	// Initialize Gmsh and create geometry
  int ierr;
  gmshInitialize(argc, argv, 0, 0, &ierr);

  int nbr_params = 4;
  double * params = (double*) malloc(nbr_params*sizeof(double));
  params[0] = 0.006; params[1] = 0.011; params[2] = 0.038; params[3] = 0.080;

  double * u; int n_nodes;
  double freq;

  double m_target = 0.050;      // 50 grams
  double freq_target = 1453.07; // F6#
  double rho = 3000;
  double v_target = design_vol(m_target, rho);
  double tol = 1; // tolerance at 1Hz

  find_param(params,3,20e-6,1e-1,tol, &u, &n_nodes, v_target, &freq, freq_target);

  if(animation) animate_in_gmsh(u, n_nodes, freq, animation_filename);

  gmshFinalize(&ierr);

  // write parameters in file
  FILE * file = fopen(results_filename, "w"); 
  fprintf(file, "optimal parameters:\nr1=%3e r2=%3e e=%3e l=%3e", params[0], params[1], params[2], params[3]);
  fclose(file);

  free(u);
  free(params);

  return 0;
}
