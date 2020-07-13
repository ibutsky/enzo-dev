/***********************************************************************
/
/  GRID CLASS (Compute and apply cold gas subgrid model source terms )
/
/  written by:  Iryna Butsky
/  date:        July, 2020
/
/  PURPOSE:  Calculates and applies cold gas subgrid model source terms
/  
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h> 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"


int grid::ComputeColdGasSourceTerms(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;


  // Some locals
  int size = 1, idx, i,j,k;
  float dvx, dvy, dvz, v2; 

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  float *drag_coef = new float[size];
  float *peak_cooling_time = new float[size];

  
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum,
      CDensNum, CVel1Num, CVel2Num, CVel3Num;


  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                   Vel3Num, TENum);
  this->IdentifyColdGasPhysicalQuantities(CDensNum, CVel1Num, CVel2Num, CVel3Num);


  // calculate cooling time at peak of cooling curve
  if (this->ComputeCoolingTimeAtSpecifiedTemperature(peak_cooling_time, CGSMMaximumCoolingTemperature, TRUE) == FAIL) {
     ENZO_FAIL("Error in grid->ComputeCoolingTimeAtSpecifiedTemperature.\n");
  }

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	idx = ELT(i,j,k);

	// First we update momentum change through cold gas drag
	if (CGSMDragModel > 0) {
	  // calculate cold gas drag coefficient
	  if (this->CalculateColdGasDragCoefficient(drag_coef) == FAIL){
		ENZO_FAIL("Error in CalculateColdGasDragCoefficient\n");
	  }
	  
	  // subtract out old kinetic energy from total energy
	  v2 = BaryonField[Vel1Num][idx] * BaryonField[Vel1Num][idx];
	  if (GridRank > 1) v2 += BaryonField[Vel2Num][idx]*BaryonField[Vel2Num][idx];
	  if (GridRank > 2) v2 += BaryonField[Vel3Num][idx]*BaryonField[Vel3Num][idx];
	  BaryonField[TENum][idx] -= 0.5 * v2;

	  if (BaryonField[CDensNum][idx] > 0) {
	    dvx = BaryonField[Vel1Num][idx] - BaryonField[CVel1Num][idx];
	    BaryonField[Vel1Num][idx]  -= dtFixed * drag_coef[idx] * dvx / BaryonField[DensNum][idx];
	    BaryonField[CVel1Num][idx] += dtFixed * drag_coef[idx] * dvx / BaryonField[CDensNum][idx];
	    if (GridRank > 1){
	      dvy = BaryonField[Vel2Num][idx] - BaryonField[CVel2Num][idx];
	      BaryonField[Vel2Num][idx]  -= dtFixed * drag_coef[idx] * dvy / BaryonField[DensNum][idx];
	      BaryonField[CVel2Num][idx] += dtFixed * drag_coef[idx] * dvy / BaryonField[CDensNum][idx];
	    }
	    if (GridRank > 2){
	      dvz = BaryonField[Vel3Num][idx] - BaryonField[CVel3Num][idx];
	      BaryonField[Vel3Num][idx]  -= dtFixed * drag_coef[idx] * dvz / BaryonField[DensNum][idx];
	      BaryonField[CVel3Num][idx] += dtFixed * drag_coef[idx] * dvz / BaryonField[CDensNum][idx];
	    }
	  }  
	  // add back the updated kinetic energy
	  v2 = BaryonField[Vel1Num][idx] * BaryonField[Vel1Num][idx];
	  if (GridRank > 1) v2 += BaryonField[Vel2Num][idx]*BaryonField[Vel2Num][idx];
	  if (GridRank > 2) v2 += BaryonField[Vel3Num][idx]*BaryonField[Vel3Num][idx];
	  BaryonField[TENum][idx] += 0.5 * v2;
	}

	// cold gas formation source terms
	if (CGSMThermalInstability > 0){
	  
	  // TODO

	  
	}
	// next: cold gas destruction source terms
	if (CGSMCloudCrushing > 0){
	  
	  // TODO
	}
  } // end triple for

  delete [] drag_coef;
  delete [] peak_cooling_time;

  return SUCCESS;  
}

