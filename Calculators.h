#ifndef _CAL_H
#define _CAL_H

#include <Mesh.h>


#include <stdio.h>      
#include <math.h>       

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
//#include "transportSimplex.h"

#define PI 3.14159265

class CCalculators
{
public:
	CCalculators();
	~CCalculators();
	
	typedef std::vector<float> nBin;
	
	bool Cal_ShapeIndex(float* cmin, float* cmax, float* SI, int n_vertex);
	std::vector<float> Make_Histogram(float* In, float bucket_size, float TEST_N);

	float Minimum(float* A, int size);
	bool Save_Ascii(float* mWrite,const char* filename, int n_vertex);
	bool Save_KWMValue(float* mWrite,const char* filename, int n_vertex);
	bool CAL_EMD(int n_BasicShape,int Test_N, std::vector<float> Basis_X, std::vector<float> SI_HIST, float* temp_EMD );
		
	bool GoldenSection(int n_BasicShape,int Test_N, std::vector<float> Basis_X, std::vector<float> SI_HIST, float* temp_EMD,int STEP);
	int OPT_EMD(float* temp_EMD, int START,std::vector<float> BASE, std::vector<float> HIST_A2, std::vector<float> HIST_A4,std::vector<float> HIST_B,int n_BasicShape);

	double Cal_SA(Mesh* mesh, float* SA);
	float Distance_TwoPoint(std::vector<float> pt1, std::vector<float> pt2);
	float Cal_TRI(float edge1,float edge2,float edge3);


};

#endif

