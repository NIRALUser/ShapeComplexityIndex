/***************************************************************************
 *   Copyright (C) 2015 by Sun Hyung Kim and Martin Styner		 					   *
 *   NeuroImage Analysis and Research Lab				  												 *
 *   Dept. of Psychiatry, University of North Carolina at Chapel Hill      *
 *   shykim@email.unc.edu	                             			               *
 *                                                                         *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <Mesh.h>
#include <Util/SurfaceUtil.h>
#include "Geodesic/Geodesic.h"
#include "Calculators.h"
#include "ComplexCLP.h"


#define BOOL bool
#define TRUE true
#define FALSE false
#define TEMPLATE_SA 74291 //SA of surf_reg_model_left.obj
//#define TEMPLATE_SA 48410.9 //SA of AVGmid_LEFT.obj, 12month IBIS data set.


int main(int argc, char *argv[])
{
	
	PARSE_ARGS;
	if(argc<2)
	{
		cout << "Cal_Complex --help" << endl;
		return false;
	}else{
		cout << "Input: " << InputSurface.c_str() << endl;	
		cout << "# of bin: " << pow(2,StepSize) + 1 << endl;
	}
	
	if(InputTemplateSurface.size() == 0 )
	{
		InputTemplateSurface = "/Users/shykim/Working_Data/Cal_Complex_Cversion/OPT_CLP_Ver/Build/template.obj";
	}
	
	CCalculators* Cal = new CCalculators;
	Mesh *MNI_GetMesh = new Mesh();
	MNI_GetMesh->openFile(InputSurface.c_str());
	
	//template surface
	Mesh *Template = new Mesh();
	Template->openFile(InputTemplateSurface.c_str());
	int tn_vertex = Template->nVertex();
	int tn_mesh = Template -> nFace();
	float *T_SA = new float[tn_vertex];
	float SUM_T_SA = 0.0f;
	float MEAN_T_SA = 0.0f;
	SUM_T_SA = Cal->Cal_SA(Template,T_SA);
	MEAN_T_SA = SUM_T_SA/tn_vertex;
	
	// NN smoothing
	SurfaceUtil::smoothing(MNI_GetMesh, 2);
	
	// principal curvature
	int n_vertex = MNI_GetMesh->nVertex();
	int n_mesh = MNI_GetMesh->nFace();
	
	float *cmin = new float[n_vertex];
	float *cmax = new float[n_vertex];
	SurfaceUtil::curvature(MNI_GetMesh, cmin, cmax, NULL, NULL);	// if principal directions are unnecessary
		
	// Surface Area
	float* SA = new float[n_vertex];
	float SUM_SA = 0.0f;
	float MEAN_SA = 0.0f;
	SUM_SA = Cal->Cal_SA(MNI_GetMesh,SA);
	MEAN_SA = SUM_SA/n_mesh;  	
	
	// calculate shape index
	float* SI = new float[n_vertex];
	Cal->Cal_ShapeIndex(cmin, cmax, SI, n_vertex);
	
	// make basic shape histogram
	int STEP = StepSize;
	int n_BasicShape =  pow(2,STEP) + 1;
  float step_size = 1/pow(2,(STEP-1));
  
  std::vector<float> Basis_X(n_BasicShape,0);
  for(int i=0; i<n_BasicShape; i++)
  {
  	Basis_X[i] = (-1) + i*step_size;
  }
  
	// set up, geodesic distance
	Geodesic geo(MNI_GetMesh);
	double dmax;
	const double *dist = geo.dist();
	std::vector<float> Testing_SI;
	
	float* EMD = new float[n_vertex];
	float* temp_EMD_Ldegree = new float[n_BasicShape];
	float* temp_EMD_Hdegree = new float[5];	
	
	//set kernel size
	if(scaleLSA == false && scaleGSA == false)
	{
		dmax = KernelSize;
		//std::cout << "dmax:" << dmax << std::endl;
	}
	if(scaleGSA == true)
	{
		dmax = KernelSize*(SUM_SA/SUM_T_SA);
		//std::cout << "dmax:"<< dmax << std::endl;
	}	
	
	for(int i=0; i<n_vertex; i++)
	//for(int i=0; i<1; i++)
	{
		if(scaleLSA == true)
		{
			dmax = KernelSize*(SA[i]/T_SA[i]);
			//std::cout << "dmax:"<< dmax << std::endl;
		}
		if(scaleGSA == true &&  scaleLSA== true)
		{
			dmax = KernelSize * (SUM_SA/SUM_T_SA) * (SA[i]/T_SA[i]);
		}			
		
		Testing_SI.clear();
		// max distance: compute geodesic distances within a specific range
		geo.perform_front_propagation(i, dmax);
		for (int j = 0; j < n_vertex; j++)
		{
			if (dist[j] <= dmax && j != i)
				{
					Testing_SI.push_back(SI[j]);
				}
		}
		
		int Test_N = Testing_SI.size()-1;
		std::vector<float> SI_HIST=Cal->Make_Histogram(&Testing_SI[0], step_size, Test_N);
		
		if (STEP<3)
		{
			Cal->CAL_EMD(n_BasicShape,Test_N, Basis_X, SI_HIST, temp_EMD_Ldegree );
			EMD[i] = 2*Cal->Minimum(temp_EMD_Ldegree, n_BasicShape);
		}else{
			Cal->GoldenSection(n_BasicShape,Test_N, Basis_X, SI_HIST, temp_EMD_Hdegree, STEP);
			EMD[i] = 2*Cal->Minimum(temp_EMD_Hdegree, 5);
		}
	}
	
	if(format == "ASCII")
	{
		Cal->Save_Ascii(EMD,output.c_str(),n_vertex);
	}
	if(format == "KWM")
	{
		Cal->Save_KWMValue(EMD,output.c_str(),n_vertex);
	} 
	
	delete [] cmin; cmin = NULL;
	delete [] cmax; cmax = NULL;
	delete [] SI; SI = NULL;
	delete [] EMD; EMD = NULL;
	delete [] temp_EMD_Ldegree; temp_EMD_Ldegree=NULL;
	delete [] temp_EMD_Hdegree; temp_EMD_Hdegree=NULL;
	delete [] SA; SA = NULL;	
	delete [] T_SA; T_SA = NULL;
	
	return true;
}

