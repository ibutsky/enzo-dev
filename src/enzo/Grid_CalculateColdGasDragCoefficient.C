/***********************************************************************
/
/  GRID CLASS (Calculates the cold gas drag terms)
/
/  written by:  Iryna Butsky
/  date:        July, 2020
/
/  PURPOSE:  Calculates cold gas drag coefficient
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


int grid::CalculateColdGasDragCoefficient(float *drag_coef){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;


  // Some locals
  int size = 1, idx, i,j,k;
  float dvx, dvy, dvz, v2, m_cold_cloudlet;
  
  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  

  float *cold_gas_density = new float[size];
  float *cold_gas_radius  = new float[size];
  
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
  float drag_units = DensityUnits / TimeUnits;
 
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	idx = ELT(i,j,k);

	// Calculate drag coefficient
	// constant drag
	if (CGSMDragModel == 1) {
	  drag_coef[idx] = CGSMDragCoefficient / drag_units;
	}

	// drag proportional to cross sectional area of cloudlets
	else if (CGSMDragModel == 2){	  
	  if (this->CalculateColdGasDensityAndRadius(cold_gas_density, cold_gas_radius) == FAIL){
	    ENZO_FAIL("Error in grid->ComputeColdGasDensityAndRadius\n");
	  }
	  m_cold_cloudlet = cold_gas_density[idx] * (4.0*PI*POW(cold_gas_radius[idx], 2)/3.0);

	  // using eq. 41 from Laibe + Price 2012b 
	  drag_coef[idx] = PI * BaryonField[DensNum][idx] * POW(cold_gas_radius[idx],2) *
	    BaryonField[CDensNum][idx] / m_cold_cloudlet; 	    
	  
	}

	else drag_coef[idx] = 0;

  } // end triple for

  delete [] cold_gas_density;
  delete [] cold_gas_radius;
  return SUCCESS;  
}

