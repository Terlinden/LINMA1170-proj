
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>

void designHalfTuningFork(double r1, double r2, double e, double l, double meshSizeFactor, char * filename) {
  /**
   * This function designs half of a tuning fork
   * r1 = inner radius (half-distance between prongs)
   * r2 = outer radius (half-width of fork)
   * e  = length of handle
   * l  = length of prongs
   * meshSizeFactor = meshSize / width of prongs
   * if `filename` is not NULL, save to file
  */
  
  int ierr;
  gmshClear(&ierr);

  double h = r2 - r1;
  double meshSize = h * meshSizeFactor;

  // Add points
  double x = 0;
  double y = 0;
  double z = 0;
  gmshModelOccAddPoint(x,y,z,meshSize,1,&ierr);
  x += h/2;
  gmshModelOccAddPoint(x,y,z,meshSize,2,&ierr);
  y += e;
  gmshModelOccAddPoint(x,y,z,meshSize,3,&ierr);
  x += r2;
  y += r2;
  gmshModelOccAddPoint(x,y,z,meshSize,4,&ierr);
  y += l;
  gmshModelOccAddPoint(x,y,z,meshSize,5,&ierr);
  x -= h;
  gmshModelOccAddPoint(x,y,z,meshSize,6,&ierr);
  y -= l;
  gmshModelOccAddPoint(x,y,z,meshSize,7,&ierr);
  x -= r1;
  y -= r1;
  gmshModelOccAddPoint(x,y,z,meshSize,8,&ierr);
  x -= h/2;
  gmshModelOccAddPoint(x,y,z,meshSize,9,&ierr);
  x += -r1 -h/2 + r2;
  y += 2*r1 - r2 + h;
  gmshModelOccAddPoint(x,y,z,meshSize,10,&ierr);
  
  // Add curves
  gmshModelOccAddLine(1,2,1,&ierr);
  gmshModelOccAddLine(2,3,2,&ierr);
  gmshModelOccAddCircleArc(3,10,4,3,&ierr);
  gmshModelOccAddLine(4,5,4,&ierr);
  gmshModelOccAddLine(5,6,5,&ierr);
  gmshModelOccAddLine(6,7,6,&ierr);
  gmshModelOccAddCircleArc(7,10,8,7,&ierr);
  gmshModelOccAddLine(8,9,8,&ierr);
  gmshModelOccAddLine(9,1,9,&ierr);

  // Add wire (closed curve)
  int curveTags[10] = {1,2,3,4,5,6,7,8,9,10};
  gmshModelOccAddWire(curveTags, 9, 1, 1, &ierr);
  
  // Add surface
  int wireTags[1] = {1};
  gmshModelOccAddPlaneSurface(wireTags, 1, 100, &ierr);

  // Sync
  gmshModelOccSynchronize(&ierr);

  // Create physical group for surface
  int surfaceTags[1] = {100};
  gmshModelAddPhysicalGroup(2, surfaceTags, 1, -1, "bulk", &ierr);

  // Create physical group for clamped curves
  int clampedCurveTags[3] = {1, 2, 9};
  gmshModelAddPhysicalGroup(1, clampedCurveTags, 3, -1, "clamped", &ierr);

  gmshModelMeshGenerate(2, &ierr);

  if(filename != NULL) gmshWrite(filename, &ierr);
}