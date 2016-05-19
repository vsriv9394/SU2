#include "../include/output_structure.hpp"

enum BCVAR  { bcMach, bcTemp, bcPres, bcDens };

void COutput::SetInriaRestart(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
	
  unsigned short nZone = geometry->GetnZone();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  bool dynamic_fem = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  //ofstream restart_file;
  string filename;
  
	unsigned long OutSol,i, npoin = geometry->GetGlobal_nPointDomain();
	int VarTyp[GmfMaxTyp];
	double bufDbl[GmfMaxTyp];
	char OutNam[1024], BasNam[1024];
	char *ptr=NULL;
	
		
  /*--- Retrieve filename from config ---*/
  
  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem){
    filename = config->GetRestart_FEMFileName();
  } else {
    filename = config->GetRestart_FlowFileName();
  }
	
	/*--- Get output name *.solb ---*/
	
	strcpy(BasNam, filename.c_str());
	ptr = strstr(BasNam,".dat");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
	ptr = strstr(BasNam,".solb");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
	sprintf(OutNam, "%s.solb", BasNam);
	
  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone);

  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(val_iZone));
  } else if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  } else if ((fem) && (config->GetWrt_Dynamic())) {
	filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  }

  /*--- Open the restart file and write the solution. ---*/
	
	OutSol = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
	
	if ( !OutSol ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to open %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}
	
  /*--- Write the restart file ---*/

	for (iVar = 0; iVar < nVar_Total; iVar++) {
		VarTyp[iVar] = GmfSca;
	}
	
	npoin = geometry->GetGlobal_nPointDomain();
	
	//printf("SET KWD : %d %d %d %d\n", OutSol, GmfSolAtVertices, npoin, nVar_Total);
	//printf("VarTyp = ");
	//for (iVar = 0; iVar < nVar_Total; iVar++) 
	//	printf("%d ",VarTyp[iVar]);
	//printf("\n");
	
	if ( !GmfSetKwd(OutSol, GmfSolAtVertices, npoin, nVar_Total, VarTyp) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to write %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}

  for (iPoint = 0; iPoint < npoin; iPoint++) {
	
    /*--- Loop over the variables and write the values to file ---*/
    for (iVar = 0; iVar < nVar_Total; iVar++) {
			bufDbl[iVar] = SU2_TYPE::GetValue(Data[iVar][iPoint]);
    }

		GmfSetLin(OutSol, GmfSolAtVertices, bufDbl);
	}
  
	if ( !GmfCloseMesh(OutSol) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Cannot close solution file %s.", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}  
	
}




/*

	Write solutions of interest : mach, dens, pres etc.

*/
void COutput::WriteInriaOutputs(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned short val_iZone) {
  
  /*--- Local variables ---*/
	
  unsigned short nZone = geometry->GetnZone();
  unsigned short Kind_Solver  = config->GetKind_Solver();
  unsigned short iVar, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iExtIter = config->GetExtIter();
  bool grid_movement = config->GetGrid_Movement();
  bool dynamic_fem = (config->GetDynamic_Analysis() == DYNAMIC);
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  //ofstream restart_file;
  string filename;
  
	unsigned long OutSol,i, npoin = geometry->GetGlobal_nPointDomain();
	int VarTyp[GmfMaxTyp];
	double bufDbl[GmfMaxTyp];
	char OutNam[1024], BasNam[1024];
	char *ptr=NULL;
	
	int NbrVar, idxVar;
	
	/* Get indices of mach, pres, etc. in the solution array */
	unsigned short *TagBc;
	TagBc = new unsigned short [100];
	
	idxVar=0;
  idxVar += nVar_Consv;
  
  if (!config->GetLow_MemoryOutput()) {
    
    if (config->GetWrt_Limiters()) {
      idxVar += nVar_Consv;
    }
    if (config->GetWrt_Residuals()) {
      idxVar += nVar_Consv;
    }
		
    if (config->GetKind_Regime() == FREESURFACE) {
      TagBc[bcDens] = idxVar;
      idxVar++;
    }
    
    if ((Kind_Solver == EULER) || (Kind_Solver == NAVIER_STOKES) || (Kind_Solver == RANS)) {
      TagBc[bcPres] = idxVar;
      TagBc[bcTemp] = idxVar+1;
      TagBc[bcMach] = idxVar+3;
      idxVar += 4;
    }
		
	}
	
  /*--- Retrieve filename from config ---*/
  
  if ((config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint())) {
    filename = config->GetRestart_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem){
    filename = config->GetRestart_FEMFileName();
  } else {
    filename = config->GetRestart_FlowFileName();
  }
	
	/*--- Get output name *.solb ---*/
	
	sprintf(OutNam, "mach.solb");
	
  /*--- Append the zone number if multizone problems ---*/
  if (nZone > 1)
    filename= config->GetMultizone_FileName(filename, val_iZone);
	
  /*--- Open the restart file and write the solution. ---*/
	
	OutSol = GmfOpenMesh(OutNam,GmfWrite,GmfDouble,nDim);
	
	if ( !OutSol ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to open %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}
	
  /*--- Write the restart file ---*/



	//for (iVar = 0; iVar < nVar_Total; iVar++) {
	//	VarTyp[iVar] = GmfSca;
	//}
	
	npoin = geometry->GetGlobal_nPointDomain();
		
	NbrVar = 1;
	VarTyp[0]  = GmfSca;
	
	if ( !GmfSetKwd(OutSol, GmfSolAtVertices, npoin, NbrVar, VarTyp) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Unable to write %s", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}
	
  for (iPoint = 0; iPoint < npoin; iPoint++) {
		
    /*--- Loop over the variables and write the values to file ---*/
    for (iVar = 0; iVar < nVar_Total; iVar++) {
	
			if ( iVar == TagBc[bcMach] ) {
				bufDbl[0] = SU2_TYPE::GetValue(Data[iVar][iPoint]);
				GmfSetLin(OutSol, GmfSolAtVertices, bufDbl);
			}
    }
		

	}
  
	if ( !GmfCloseMesh(OutSol) ) {
	  printf("\n\n   !!! Error !!!\n" );
    printf("Cannot close solution file %s.", OutNam);
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
	}  
	
	delete [] TagBc;
	
}



