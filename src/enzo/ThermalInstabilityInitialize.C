/***********************************************************************
/
/   Thermal Instability Test Problem
/
/   written by: Cameron Hummels, Iryna Butsky
/   date:       March 2018
/
/   PURPOSE: Investigate how thermal instability operates when different
/            physics are present
/
/   RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
 
int ThermalInstabilityInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
          TopGridData &MetaData)
{
   if(debug){
      printf("Entering ThermalInstabilityInitialize\n");
      fflush(stdout);
      }

   // Field Names
   char *DensName = "Density";
   char *TEName   = "TotalEnergy";
   char *GEName   = "GasEnergy";
   char *Vel1Name = "x-velocity";
   char *Vel2Name = "y-velocity";
   char *Vel3Name = "z-velocity";
   char *BxName   = "Bx";
   char *ByName = "By";
   char *BzName = "Bz";
   char *PhiName = "Phi";
   char *Phi_pName = "Phip";
   char *CRName = "CREnergyDensity";
   char *ElectronName = "Electron_Density";
   char *HIName = "HI_Density";
   char *HIIName = "HII_Density";
   char *HeIName = "HeI_Density";
   char *HeIIName = "HeII_Density";
   char *HeIIIName = "HeIII_Density";
   char *HMName = "HM_Density";
   char *H2IName = "H2I_Density";
   char *H2IIName = "H2II_Density";
   char *DIName = "DI_Density";
   char *DIIName = "DII_Density";
   char *HDIName = "HDI_Density";
   char *MetalName = "Metal_Density";
   char *ExternalAccelxName = "ExternalAcceleration_x";
   char *ExternalAccelyName = "ExternalAcceleration_y";
   char *ExternalAccelzName = "ExternalAcceleration_z";


   /* parameter declarations */
 
   FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

   // Local variable declarations
   char line[MAX_LINE_LENGTH];
   int   dim, ret, level;

   // Initialize parameters to default values
   float TIMeanDensity = 1e-28;
   float TIPerturbationAmplitude = 0.01;
   float TIMeanTemperature = 1e6;
   int TIHaloProfile = 0;
   int TIPerturbationType = 1;
   // Where InverseBeta = P_magnetic / P_gas to easily allow for B = 0
   int TIMagneticFieldUseConstantBeta = 1;
   int TICosmicRayUseConstantEta = 1;
   float TIMagneticFieldInverseBeta = 0;
   float TIMagneticFieldDirection[MAX_DIMENSION];
   // P_cr / P_gas
   float TICosmicRayPressureRatio = 0;
   for (int i = 0; i < MAX_DIMENSION; i++)
     TIMagneticFieldDirection[i] = 0; 

   // Read problem specific parameters. 
   while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
      ret = 0;

      ret += sscanf(line, "TIMeanDensity = %"FSYM, &TIMeanDensity);
      ret += sscanf(line, "TIPerturbationAmplitude = %"FSYM, &TIPerturbationAmplitude);
      ret += sscanf(line, "TIMeanTemperature = %"FSYM, &TIMeanTemperature);
      ret += sscanf(line, "TIHaloProfile = %"ISYM, &TIHaloProfile);
      ret += sscanf(line, "TIPerturbationType = %"ISYM, &TIPerturbationType);
      ret += sscanf(line, "TestProblemUseMetallicityField = %"ISYM, &TestProblemData.UseMetallicityField);
      ret += sscanf(line, "TIMagneticFieldDirection = %"FSYM" %"FSYM" %"FSYM,
                    &TIMagneticFieldDirection[0], &TIMagneticFieldDirection[1], &TIMagneticFieldDirection[2]);
      ret += sscanf(line, "TIMagneticFieldUseConstantBeta = %"ISYM, &TIMagneticFieldUseConstantBeta);
      ret += sscanf(line, "TICosmicRayUseConstantEta = %"ISYM, &TICosmicRayUseConstantEta);
      ret += sscanf(line, "TIMagneticFieldInverseBeta = %"FSYM, &TIMagneticFieldInverseBeta);
      ret += sscanf(line, "TICosmicRayPressureRatio = %"FSYM, &TICosmicRayPressureRatio);

      // Issue a warning if the line is suspicious 
      if (ret == 0 && strstr(line, "=") && strstr(line, "TI") && line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
         fprintf(stderr, "*** warning: the following parameter line was not interpreted:\n%s\n", line);
    } // end input from parameter file

    TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)

    /* set the boundary conditions */

    /*    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
        MetaData.LeftFaceBoundaryCondition[dim]  = LeftFaceBoundaryCondition[dim];
        MetaData.RightFaceBoundaryCondition[dim] = RightFaceBoundaryCondition[dim];
    }*/

    /* Initialize a uniform grid first */

   float uniform_velocity[3] = {0.0, 0.0, 0.0};
   float uniform_density = 1.0;
   float uniform_total_energy = 1.0;
   float uniform_B_field[3] = {0.0, 0.0, 0.0};

   if (TopGrid.GridData->InitializeUniformGrid(uniform_density,
                        uniform_total_energy,
                        uniform_total_energy,
                        uniform_velocity,
                        uniform_B_field) == FAIL) {
                        ENZO_FAIL("Error in InitializeUniformGrid.");
                        }
    
    /* set up grid fields from turbulent ICs*/
   if (TopGrid.GridData->ThermalInstabilityInitializeGrid(TIMeanDensity, TIMeanTemperature,
                                                           TIPerturbationAmplitude,
							   TIPerturbationType, TIHaloProfile, 
							   TIMagneticFieldDirection, 
							   TIMagneticFieldUseConstantBeta,
 							   TICosmicRayUseConstantEta, 
							   TIMagneticFieldInverseBeta, 
							   TICosmicRayPressureRatio) 
							    == FAIL) {
        ENZO_FAIL("Error in ThermalInstabilityInitializeGrid.");
    }

    /* set up field names and units */

    int i = 0;
    int j;
    DataLabel[i++] = DensName;
    DataLabel[i++] = TEName;
    if (DualEnergyFormalism)
        DataLabel[i++] = GEName;
    DataLabel[i++] = Vel1Name;
    DataLabel[i++] = Vel2Name;
    DataLabel[i++] = Vel3Name;
    if (HydroMethod == MHD_RK) {
      DataLabel[i++] =  BxName;
      DataLabel[i++] =  ByName;
      DataLabel[i++] =  BzName;
      DataLabel[i++] =  PhiName;
    }
    if (UsePoissonDivergenceCleaning){
      DataLabel[i++] = Phi_pName;
    }
    if (CRModel)
      DataLabel[i++] = CRName;

    if (TestProblemData.MultiSpecies) {
      DataLabel[i++] = ElectronName;
      DataLabel[i++] = HIName;
      DataLabel[i++] = HIIName;
      DataLabel[i++] = HeIName;
      DataLabel[i++] = HeIIName;
      DataLabel[i++] = HeIIIName;

      if (TestProblemData.MultiSpecies > 1) {
         DataLabel[i++] = HMName;
         DataLabel[i++] = H2IName;
         DataLabel[i++] = H2IIName;
       }

      if (TestProblemData.MultiSpecies > 2) {
         DataLabel[i++] = DIName;
         DataLabel[i++] = DIIName;
         DataLabel[i++] = HDIName;
       }
    }
    if (TestProblemData.UseMetallicityField)
        DataLabel[i++] = MetalName;

   for(j=0; j < i; j++)
      DataUnits[j] = NULL;

    /* Write parameters to parameter output file */

    if (MyProcessorNumber == ROOT_PROCESSOR) {
        fprintf(Outfptr, "TIMeanDensity  = %"FSYM"\n", TIMeanDensity);
        fprintf(Outfptr, "TIPerturbationAmplitude = %"FSYM"\n", TIPerturbationAmplitude);
        fprintf(Outfptr, "TIMeanTemperature  = %"FSYM"\n", TIMeanTemperature);
	fprintf(Outfptr, "TIHaloProfile = %"ISYM"\n", TIHaloProfile);
    }

    return SUCCESS;

}
