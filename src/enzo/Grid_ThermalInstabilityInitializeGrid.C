/***********************************************************************
/
/   GRID CLASS (INITIALIZE THE GRID FOR THE THERMAL INSTABILITY TEST)
/
/   written by: Iryna Butsky and Cameron Hummels
/   date:         June 2018
/   modified1:   
/
/   PURPOSE: Sets up the grid for the ThermalInstability problem type.
/            Your run directory must include the file 'perturbation.in', 
/            which can be generated using the file 'white_noise_generator.py', 
/            which should be in run/Hydro/Hydro-3D/ThermalInstability
/
/   RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::ThermalInstabilityInitializeGrid(float TIMeanDensity, float TIMeanTemperature, 
					   float TIPerturbationAmplitude,
					   int TIPerturbationType, int TIHaloProfile, 
					   float *TIMagneticFieldDirection, 
					   int TIMagneticFieldUseConstantBeta,
					   int TICosmicRayUseConstantEta, 
					   float TIMagneticFieldInverseBeta, 
					   float TICosmicRayPressureRatio)
{
   if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

   if(debug){
      printf("Entering ThermalInstabilityInitializeGrid\n");
      fflush(stdout);
      }
 
   // Figure out grid quantities and how to access fields 
   int size = 1;

   for (int dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

   int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num, PhiNum, CRNum;
   int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, MetalNum;

   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
      }

   if (CRModel) {
     if ((CRNum = FindField(CRDensity, FieldType, NumberOfBaryonFields)) < 0)
       ENZO_FAIL("Cannot Find Cosmic Rays");
   }

   int MetallicityField = FALSE;

   if (MultiSpecies) {
      if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
         ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
         }
      }

   if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
         != -1)
      MetallicityField = TRUE;
   else
      MetalNum = 0;

   // Get the units
   float AccelerationUnits, DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits, 
     MassUnits;

   
   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	    &TimeUnits, &VelocityUnits, &MassUnits, Time);

   AccelerationUnits = LengthUnits / TimeUnits / TimeUnits;

   printf("Getting ready to set up the thermal instability box.\n");

   // Set up the white noise field.
   int pert_size_x;
   int pert_size_y;
   int pert_size_z;

   float*** perturbation_field;

   
   FILE* inf;
   inf = fopen("perturbation.in", "r");

   if (inf == NULL) {
      printf("Could not open file 'perturbation.in', which is needed for the pressure perturbation.\n");
      printf("If it is not here, try running the file 'white_noise_generator.py' in run/ dir.\n");
      exit(1);
    }
      
      // Read in comments and the first two lines (the dimension and the grid size)
      int lines_read = 0;
      size_t line_length = 80;
      
      int pert_dim;
      float pert_std;

      char* line = NULL;

      char* pert_dim_x_s = new char[line_length];
      char* pert_dim_y_s = new char[line_length];
      char* pert_dim_z_s = new char[line_length];
      char* pert_std_s   = new char[line_length];
      while (lines_read < 3) {
         getline(&line, &line_length, inf);
         
         if (line[0] != '#') {
            if (lines_read == 0) {
               sscanf(line, "%i", &pert_dim);
               lines_read ++;
	    }

            else if (lines_read == 1){
               sscanf(line, "%s %s %s", pert_dim_x_s, pert_dim_y_s, pert_dim_z_s);
               lines_read ++;
	    }
	    else {
	      sscanf(line, "%s", pert_std_s);
	      lines_read ++;
	    }
	 } 
      }
      
     
      pert_size_x = atoi(pert_dim_x_s);
      pert_size_y = atoi(pert_dim_y_s);
      pert_size_z = atoi(pert_dim_z_s);
      pert_std    = atof(pert_std_s);
      printf("%f\n", pert_std);
      printf("Reading in a %i dim perturbation grid of size %i x %i x %i\n", pert_dim, pert_size_x, pert_size_y, pert_size_z);

      // Create an array to hold the turbulence
      perturbation_field = new float**[pert_size_z];

      for (int j = 0; j < pert_size_z; j++) {
	perturbation_field[j] = new float*[pert_size_y];

         for (int i = 0; i < pert_size_y; i++){
            perturbation_field[j][i] = new float[pert_size_x];
          }
       }

      if (GridRank < 3)
	pert_size_z = 1;

      for (int k = 0; k < pert_size_z; k++)
         for (int j = 0; j < pert_size_y; j++)
            for (int i = 0; i < pert_size_x; i++)
               perturbation_field[k][j][i] = 0.0;
       

      // Read in the turbulence.
      char* x_s = new char[line_length];
      char* y_s = new char[line_length];
      char* z_s = new char[line_length];

      char* perturbation_s = new char[line_length];

      while (getline(&line, &line_length, inf) != -1) {
         int x, y, z;

         sscanf(line, "%s %s %s %s", x_s, y_s, z_s, perturbation_s);
         x = atoi(x_s);
         y = atoi(y_s);
         z = atoi(z_s);

         perturbation_field[z][y][x] = atof(perturbation_s);
         }
      
      printf("Done reading in white noise.\n");

      for (int k = 0; k < pert_size_z; k++)
         for (int j = 0; j < pert_size_y; j++)
            for (int i = 0; i < pert_size_x; i++) 
               perturbation_field[k][j][i] *= (TIPerturbationAmplitude / pert_std);


      // Free some memory.
      delete [] x_s;
      delete [] y_s;
      delete [] z_s;

      delete [] pert_dim_x_s;
      delete [] pert_dim_y_s;
      delete [] pert_dim_z_s;

      printf("Finished reading in perturbation.\n");

   // Set grid values.
   printf("Setting grid values. This might take a minute.\n");

   int cell_index;
   float x, y, z, zscale, Bx, By, Bz, B2, halo_scale, total_pressure_factor;
   float NormalizedPerturbation, GasDensity, GasTemperature, GasPressure, dGasPressure;

   // temporary
   float RadiativeCoolingPowerLawIndex = -0.5;
   
   // define in code units
   float T0 = TIMeanTemperature / TemperatureUnits;
   float rho0 = TIMeanDensity / DensityUnits;

   float g0 = ExternalGravityConstant / AccelerationUnits;

   float a = ExternalGravitySofteningRadius * kpc_cm / LengthUnits;
   float H = kboltz * TIMeanTemperature / (Mu * mh * ExternalGravityConstant) / LengthUnits;
   float zc = ExternalGravityPosition[2]; // note: check units

   printf("a = %e, H = %e\n", a, H);
   for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
         for (int i = 0; i < GridDimension[0]; i++) {
            cell_index = k * (GridDimension[1] * GridDimension[0]) + j * GridDimension[0] + i;

            // x, y, z, radius in code units
	    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	    y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	    z = y; 
	    if (GridRank > 2)
	      z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	    zscale = (z - zc) / a;

            int pert_index[3];

            pert_index[0] = floor((float)pert_size_x * (x - DomainLeftEdge[0]) / (DomainRightEdge[0] - DomainLeftEdge[0]));
            pert_index[1] = floor((float)pert_size_y * (y - DomainLeftEdge[1]) / (DomainRightEdge[1] - DomainLeftEdge[1]));
	    if (GridRank > 2)
	      pert_index[2] = floor((float)pert_size_z * (z - DomainLeftEdge[2]) / (DomainRightEdge[2] - DomainLeftEdge[2]));
	    else 
	      pert_index[2] = 0;

            if (pert_index[0] < 0)
               pert_index[0] = 0;
            if (pert_index[1] < 0)
               pert_index[1] = 0;
            if (pert_index[2] < 0)
               pert_index[2] = 0;


	    if (pert_index[0] > pert_size_x - 1)
		pert_index[0] = pert_size_x - 1;
            if (pert_index[1] > pert_size_y - 1)
               pert_index[1] = pert_size_y - 1;
	    if (pert_index[2] > pert_size_z - 1)
		pert_index[2] = pert_size_z - 1;
	    NormalizedPerturbation = perturbation_field[pert_index[2]][pert_index[1]][pert_index[0]];

	    total_pressure_factor = 1.0;
	    
	    if (TIMagneticFieldUseConstantBeta)
	      total_pressure_factor += TIMagneticFieldInverseBeta;
	    if (TICosmicRayUseConstantEta)
	      total_pressure_factor += TICosmicRayPressureRatio;
	    halo_scale = a * (sqrt(1.0 + zscale*zscale)-1.0) / H / total_pressure_factor;

	    // First, solve for the unperturbed gas density and temperature for different profiles
	    // Uniform box, for testing
	    if (TIHaloProfile == 0) {
	      GasDensity     = rho0;
	      GasTemperature = T0; 
	    }
	    // Isothermal profile
	    else if (TIHaloProfile == 1) {
	      GasDensity     = rho0 * exp(-halo_scale);
	      GasTemperature = T0;
	    }
	    // Isentropic profile
	    else if (TIHaloProfile == 2) {
	      GasTemperature = T0 * (1.0 - ((Gamma - 1.0)/Gamma)*halo_scale);
	      GasDensity     = rho0 *  pow((1.0 - ((Gamma - 1.0)/Gamma)*halo_scale), (Gamma - 1.0));
	    }

	    // Iso t-cool profile
	    else if (TIHaloProfile == 3) {
	      GasTemperature = T0 * (1.0 - halo_scale / (2.0 - RadiativeCoolingPowerLawIndex));
	      GasDensity = rho0 * pow((1.0 - halo_scale / (2.0 - RadiativeCoolingPowerLawIndex)), 1.0 - RadiativeCoolingPowerLawIndex);
	    }
	    // note: this is gas pressure / mass
	    GasPressure = GasTemperature / Mu; 

	    // constant thermal pressure 
	    if (TIPerturbationType == 1){
	      BaryonField[DensNum][cell_index] = GasDensity * (1.0 + NormalizedPerturbation);
	      BaryonField[TENum][cell_index]   = GasPressure / (Gamma - 1.0) / (1.0 + NormalizedPerturbation);
	      if (DualEnergyFormalism)
		BaryonField[GENum][cell_index] = GasPressure / (Gamma - 1.0) / (1.0 + NormalizedPerturbation);

	    }
	    // constant gas density and constant gas + cr pressure (see below)
	    else if (TIPerturbationType == 2){
              BaryonField[DensNum][cell_index] = GasDensity;
	      dGasPressure = GasPressure * NormalizedPerturbation;
              BaryonField[TENum][cell_index]   = (GasPressure + dGasPressure) / (Gamma - 1.0);
              if (DualEnergyFormalism)
		BaryonField[GENum][cell_index] = (GasPressure + dGasPressure) / (Gamma - 1.0);
	    }

            // Set the velocity to zero 
            BaryonField[Vel1Num][cell_index] = 0.0;
            BaryonField[Vel2Num][cell_index] = 0.0;
            BaryonField[Vel3Num][cell_index] = 0.0;       
	    
	    
	    if (HydroMethod == MHD_RK){
	      if (TIMagneticFieldUseConstantBeta) 
		B2 = 2.0 * TIMagneticFieldInverseBeta * GasPressure * GasDensity; 
	      else
		B2 = 2.0 * TIMagneticFieldInverseBeta * T0 / Mu * rho0;
	      
	      Bx = sqrt(B2) * TIMagneticFieldDirection[0];
	      By = sqrt(B2) * TIMagneticFieldDirection[1];
	      Bz = sqrt(B2) * TIMagneticFieldDirection[2];

	      BaryonField[B1Num][cell_index] = Bx; 
	      BaryonField[B2Num][cell_index] = By; 
	      BaryonField[B3Num][cell_index] = Bz; 
	      BaryonField[PhiNum][cell_index] = 0.0; 

	      BaryonField[TENum][cell_index] += B2/2.0/GasDensity;
	    }

	    if (CRModel){
	      if (TICosmicRayUseConstantEta)
		BaryonField[CRNum][cell_index]  = TICosmicRayPressureRatio * GasPressure * GasDensity / (CRgamma - 1.0);	      
	      else
		BaryonField[CRNum][cell_index]  = TICosmicRayPressureRatio * T0 / Mu * rho0 / (CRgamma - 1.0);
	      if (TIPerturbationType == 2)
		BaryonField[CRNum][cell_index] -= dGasPressure / (CRgamma - 1.0) * GasDensity;
	    }

	    // Set up the chemistry.
            if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE)
	      BaryonField[MetalNum][cell_index] = BaryonField[DensNum][cell_index]*TestProblemData.MetallicityField_Fraction;
	      
         } //End loop over i
      } //End loop over j
   } //End loop over k

   printf("Finished setting grid values.\n");

   // Free the perturbation array
   for (int j = 0; j < pert_size_z; j++) {
      for (int i = 0; i < pert_size_y; i++) {
         delete[] perturbation_field[j][i];
      }
      delete[] perturbation_field[j];
   }
   delete[] perturbation_field;


   printf("All finished initializing this grid.\n");

   // Done with initialization
   if(debug){
      printf("Exiting ThermalInstabilityInitialize\n");
      fflush(stdout);
      }

   return SUCCESS;
}

