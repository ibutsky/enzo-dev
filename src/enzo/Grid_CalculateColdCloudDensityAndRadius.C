/***********************************************************************
/
/  GRID CLASS (Calculate the radius of unresolved cold cloudlets)
/
/  written by:  Iryna Butsky
/  date:        July, 2020
/
/  PURPOSE:  Calculates the radius of cold gas clouds assuming unresolved
/            cold gas is in pressure equilibrium with the hot gas phase
/            so and that r_cloud ~ min(c_s * t_cool)
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
#include "hydro_rk/EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int grid::CalculateColdCloudDensityAndRadius(float *cold_cloud_density, float *cold_cloud_radius){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;


  // Some locals
  int size = 1, idx, i,j,k;
  float v2, rho_cold, T_hot, p_hot, eint_hot, h, cs, dpdrho, dpde;

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum,
      CDensNum, CVel1Num, CVel2Num, CVel3Num;


  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                   Vel3Num, TENum);
  this->IdentifyColdGasPhysicalQuantities(CDensNum, CVel1Num, CVel2Num, CVel3Num);


  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);
  
  float *min_cooling_time = new float[size];
  
  if (this->ComputeCoolingTimeAtSpecifiedTemperature(min_cooling_time, CGSMMaximumCoolingTemperature, TRUE) == FAIL) {
      ENZO_FAIL("Error in grid->ComputeCoolingTimeAtSpecifiedTemperature.\n");
  }
  
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	idx = ELT(i,j,k);


	// cloud size ~ sound_speed * cooling_time
	// calculating sound speed in cold gas cloud assuming it's in pressure
	// equilibrium with hot phase

	if(DualEnergyFormalism)
	  eint_hot = BaryonField[GENum][idx];
	else{
	  v2 = BaryonField[Vel1Num][idx]*BaryonField[Vel2Num][idx];
	  if(GridRank > 1) v2 += BaryonField[Vel2Num][idx]*BaryonField[Vel2Num][idx];
	  if(GridRank > 2) v2 += BaryonField[Vel3Num][idx]*BaryonField[Vel3Num][idx];
	  eint_hot = BaryonField[TENum][idx] - 0.5*v2;
	}
	
	T_hot = eint_hot * Mu * (Gamma - 1.0) * TemperatureUnits;
	p_hot = eint_hot * (Gamma - 1.0) * BaryonField[DensNum][idx];
	// assuming p_cold = p_hot and ideal gas law
	// rho_cold is the true, physical density of the cold cloudlets
	// not to be confused with the volume-density of the cold gas subgrid fluid
	rho_cold  = BaryonField[DensNum][idx] * T_hot / CGSMCharacteristicTemperature;
	
	// calculate sound speed given 
	EOS(p_hot, rho_cold, eint_hot, h, cs, dpdrho, dpde, EOSType, 1);
	cold_cloud_density[idx] = rho_cold;
	cold_cloud_radius[idx] = cs * min_cooling_time[idx];
	
  } // end triple for

  delete [] min_cooling_time;
  return SUCCESS;  
  }

