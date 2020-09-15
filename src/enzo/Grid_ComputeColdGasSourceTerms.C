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
             float *VelocityUnits, FLOAT Time);

int grid::ComputeColdGasSourceTerms(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  // Some locals
  int size = 1, idx, i,j,k;
  float dvx, dvy, dvz, v2, drho, HotGasTemperature,
    dE_radiative_cooling, rho, gasenergy, coldgasenergy, dE_hot_to_cold;

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  float *dx = new float[3];

  dx[0] = CellWidth[0][0];
  dx[1] = (GridRank > 1) ? CellWidth[1][0] : 1.0;
  dx[2] = (GridRank > 2) ? CellWidth[2][0] : 1.0;

  float *drag_coef = new float[size];
  float *peak_cooling_time = new float[size];
  float *cooling_time     = new float[size];
  float *cold_gas_radius  = new float[size];
  float *cold_gas_density = new float[size];
  
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum,
      CDensNum, CVel1Num, CVel2Num, CVel3Num;
  
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                   Vel3Num, TENum);
  this->IdentifyColdGasPhysicalQuantities(CDensNum, CVel1Num, CVel2Num, CVel3Num);
 
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);

  if (CGSMDragModel > 0){
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	  idx = ELT(i,j,k);

	  // First we update momentum change through cold gas drag
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
  
      } // end for
  }
  
  // First apply momentum transfer between cold and hot gas, then energy transfer
  // Should we switch the order?

  if (CGSMThermalInstability > 0) {
    // calculate cooling time at peak of cooling curve
    if (this->ComputeCoolingTimeAtSpecifiedTemperature(peak_cooling_time, CGSMMaximumCoolingTemperature, TRUE) == FAIL) {
      ENZO_FAIL("Error in grid->ComputeCoolingTimeAtSpecifiedTemperature.\n");
    }
   if (this->ComputeCoolingTime(cooling_time, TRUE, FALSE) == FAIL) {
      ENZO_FAIL("Error in grid->ComputeCoolingTimeAtSpecifiedTemperature.\n");
    }

    if (this->CalculateColdGasDensityAndRadius(cold_gas_density, cold_gas_radius) == FAIL){
      ENZO_FAIL("Error in grid->ComputeColdGasDensityAndRadius\n");
    }
    

    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
        for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
            for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
                idx = ELT(i,j,k);

                if (DualEnergyFormalism)
                    gasenergy = BaryonField[GENum][idx];
                else {
                    v2 = BaryonField[Vel1Num][idx] * BaryonField[Vel1Num][idx];
                    if (GridRank > 1) v2 += BaryonField[Vel2Num][idx]*BaryonField[Vel2Num][idx];
                    if (GridRank > 2) v2 += BaryonField[Vel3Num][idx]*BaryonField[Vel3Num][idx];
                    gasenergy = BaryonField[TENum][idx] - 0.5*v2;
                }
                HotGasTemperature = gasenergy*pow(VelocityUnits, 2)*Mu*mh*(Gamma-1) / kboltz;
                rho = BaryonField[DensNum][idx];
                dE_radiative_cooling = rho * gasenergy * dtFixed / cooling_time[idx];

                // If gas is thermally unstable, transfer mass from hot gas to cold gas in lieu of cooling
                if ((CGSMThermalInstability > 0) & (HotGasTemperature < 1e6) & (HotGasTemperature > 1e4) & (dE_radiative_cooling < 0)){//} & (cold_gas_radius[idx] < dx[0])){
                    coldgasenergy = CGSMCharacteristicTemperature * kboltz / (Mu*mh*(Gamma-1.0)) / pow(VelocityUnits, 2);
                    // mass to be transferred
                    drho = -dE_radiative_cooling / coldgasenergy;
                    // updated hot gas energy
                    gasenergy = (rho*gasenergy + dE_radiative_cooling) / (rho - drho);
                    
                    // transfer mass between hot phase and cold phase
                    BaryonField[DensNum][idx] -= drho;
                    BaryonField[CDensNum][idx] += drho;
                    
                }
                // Else, do normal cooling
                else {
                    gasenergy += dE_radiative_cooling / rho;
                }
                
                // in both cases, need to update the internal and total energies
                if (DualEnergyFormalism)
                    BaryonField[GENum][idx] = gasenergy;
                // update the total energy with the new gas energy
                v2 = BaryonField[Vel1Num][idx] * BaryonField[Vel1Num][idx];
                if (GridRank > 1) v2 += BaryonField[Vel2Num][idx]*BaryonField[Vel2Num][idx];
                if (GridRank > 2) v2 += BaryonField[Vel3Num][idx]*BaryonField[Vel3Num][idx];
                BaryonField[TENum][idx] = gasenergy + 0.5*v2;
        
                
          } // end triple for
    } // end if CGSMThermalInstability
        
	// next: cold gas destruction source terms
	if (CGSMCloudCrushing > 0){
	  
	  // TODO
	}
    
  delete [] dx;
  delete [] drag_coef;
  delete [] peak_cooling_time;
  delete [] cooling_time;			       
  delete [] cold_gas_radius;
  delete [] cold_gas_density;
  return SUCCESS;  
}

