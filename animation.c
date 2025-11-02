#include <stdlib.h>
#include <stdio.h>
#include <gmshc.h>
#include <math.h>
#include <unistd.h>
#include "animation.h"

#define FPS 24
#define ANIMATION_PERIOD 1
#define VIDEO_LENGTH 4
#define PI 3.14

void animate_in_gmsh(double* SOL, int n_nodes, int frequency, char * name){

  char * input_folder    = "frames";
  char output_filename[100];
  sprintf(output_filename, "animations/%s", name);

  // Remove all of the pictures before running

  char delete_cmd[100];
  sprintf(delete_cmd, "rm %s/*", input_folder);

  system(delete_cmd);

  int ierr;
  double time;
  char legend[100];

  int num_frames = VIDEO_LENGTH * FPS;

  size_t *nodeTags=malloc(n_nodes*sizeof(size_t));
  
  gmshFltkInitialize(&ierr);

  for(int frame = 0; frame < num_frames; frame++) {

    double* sol_3D = malloc(3*n_nodes*sizeof(double)); // in gmsh, vectors are 3D
    for (int i=0; i<n_nodes; i++){
        nodeTags[i] = i+1;
        sol_3D[3*i+0] = SOL[2*i+0] * cos(2*PI*frame/(FPS*ANIMATION_PERIOD));
        sol_3D[3*i+1] = SOL[2*i+1] * cos(2*PI*frame/(FPS*ANIMATION_PERIOD));
        sol_3D[3*i+2] = 0.;
    }

    char *model_name;
    gmshModelGetCurrent(&model_name, &ierr);
    
    time = (double) frame / (FPS*ANIMATION_PERIOD*frequency);
    sprintf(legend, "displacement at time: %.2es.", time);
    int view_tag = gmshViewAdd(legend, -1, &ierr);
    gmshViewAddHomogeneousModelData(view_tag, 0, model_name,"NodeData", nodeTags, n_nodes, sol_3D, 3*n_nodes, 0., 3, -1, &ierr);
    gmshViewOptionSetNumber(view_tag, "VectorType", 5, &ierr);
    gmshViewOptionSetNumber(view_tag, "DisplacementFactor", 0.1, &ierr);
  
    char filename[100];
    sprintf(filename, "%s/frame%d.png", input_folder, frame);

    gmshWrite(filename, &ierr);

    gmshViewRemove(view_tag, &ierr);
  }

  gmshFltkFinalize(&ierr);

  // create the video
  char command[256];
  snprintf(command, sizeof(command),
    "ffmpeg -framerate %d -i %s/frame%%d.png -c:v libx264 -r 30 -pix_fmt yuv420p %s",
    FPS, input_folder, output_filename);

  system(command);
}