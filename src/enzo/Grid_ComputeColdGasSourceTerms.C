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

int GetUnits(float *DensityUnits, float *LengthUnits,
            float *TemperatureUnits, float *TimeUnits,
            float *VelocityUnits, double *MassUnits, FLOAT Time);


int grid::ComputeColdGasSourceTerms(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;


  // Some locals
  int size = 1, idx, i,j,k;
  float drag_coef, dvx, dvy, dvz; 

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum,
      CDensNum, CVel1Num, CVel2Num, CVel3Num;


  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                   Vel3Num, TENum);
  this->IdentifyColdGasPhysicalQuantities(CDensNum, CVel1Num, CVel2Num, CVel3Num);

 // Some locals
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0;
  double MassUnits = 1.0;

  // Get system of units
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  float drag_units = MassUnits / LengthUnits / LengthUnits / TimeUnits / TimeUnits;
  drag_coef = CGSMDragCoefficient / drag_units; 
 
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	idx = ELT(i,j,k);

	// subtract out old kinetic energy from total energy
	//	BaryonField[TENum][idx] -= 0.5 * (BaryonField[Vel1Num][idx]*BaryonField[Vel1Num][idx] + 
	//					  ((GridRank > 1) ? BaryonField[Vel2Num][idx]*BaryonField[Vel2Num][idx] : 0) +
	//					  ((GridRank > 2) ? BaryonField[Vel3Num][idx]*BaryonField[Vel3Num][idx] : 0));

	if (BaryonField[CDensNum][idx] > 0) {
	  dvx = BaryonField[Vel1Num][idx] - BaryonField[CVel1Num][idx];
	  BaryonField[Vel1Num][idx] -= dtFixed * drag_coef * dvx / BaryonField[DensNum][idx];
	  BaryonField[CVel1Num][idx] += dtFixed * drag_coef * dvx / BaryonField[CDensNum][idx];

	  if (GridRank > 1){
	    dvy = BaryonField[Vel2Num][idx] - BaryonField[CVel2Num][idx];
	    BaryonField[Vel2Num][idx] -= dtFixed * drag_coef * dvy / BaryonField[DensNum][idx];
	    BaryonField[CVel2Num][idx] += dtFixed * drag_coef * dvy / BaryonField[CDensNum][idx];
	  }
	  if (GridRank > 2){
	    dvz = BaryonField[Vel3Num][idx] - BaryonField[CVel3Num][idx];
	    BaryonField[Vel3Num][idx] -= dtFixed * drag_coef * dvz / BaryonField[DensNum][idx];
	    BaryonField[CVel3Num][idx] += dtFixed * drag_coef * dvz / BaryonField[CDensNum][idx];
	  }
	}  
	// add back the updated kinetic energy
	//        BaryonField[TENum][idx] += 0.5 * (BaryonField[Vel1Num][idx]*BaryonField[Vel1Num][idx] +
	//                          (GridRank > 1) ? BaryonField[Vel2Num][idx]*BaryonField[Vel2Num][idx] : 0 +
	//                          (GridRank > 2) ? BaryonField[Vel3Num][idx]*BaryonField[Vel3Num][idx] : 0);

  } // end triple for

  return SUCCESS;  
}

