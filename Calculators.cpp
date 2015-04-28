#include "Calculators.h"
#include "transportSimplex.h"

using namespace t_simplex;
using namespace std;


CCalculators::CCalculators()
{
	
	
}

CCalculators::~CCalculators()
{


}

double DISTANCE(float A,float B) 
{
	//Formula for approximating the distance between a factory and a warehouse;
	return std::abs(A-B);
};



bool CCalculators::Cal_ShapeIndex(float* cmin, float* cmax, float* SI, int n_vertex)
{
	
	for(int i=0; i<n_vertex; i++)
	{
		SI[i] = (2/PI)*atan( (cmax[i]+cmin[i])/(cmax[i]-cmin[i]) );
		
	}
	return true;
		
}

std::vector<float> CCalculators::Make_Histogram(float* In, float bucket_size, float TEST_N)
{
	std::vector<float> out_hist;
	
	//initial shift, range 0 to 2
	for(int i=0; i<TEST_N; i++){
		In[i] = In[i] + 1;
	}
	
	int number_of_buckets = (int)ceil(1 / bucket_size)*2; 
	std::vector<int> histogram(number_of_buckets+1);
	
	for(int i=0; i<TEST_N; i++){
	    int bucket = (int)(In[i] / bucket_size);
	    histogram[bucket] += 1;
	}
	
	for (int i=0; i<number_of_buckets+1; i++)
	{
		out_hist.push_back(histogram[i]);
	}
		
	return out_hist;	
}


float CCalculators::Minimum(float* A, int size)
{ 
	float temp = 10.0f;
	for(int i=0;i<size; i++){
		if(A[i]<temp){
			temp = A[i];
		}
	}	
		
	return temp;
}
	
bool CCalculators::Save_Ascii(float* mWrite,const char* filename, int n_vertex)
{
	FILE* writefile;
  
	writefile = fopen( filename, "w" );
	
	for(int i=0; i<n_vertex; i++){
		fprintf( writefile, "%f\n", mWrite[i]);
	}
			
	fclose(writefile);
	return true;
}		

bool CCalculators::Save_KWMValue(float* mWrite,const char* filename, int n_vertex)
{
	FILE* writefile;
  
	writefile = fopen( filename, "w" );
	
	fprintf(writefile, "NUMBER_OF_POINTS=");
	fprintf(writefile, "%i \n", n_vertex);
	fprintf(writefile, "DIMENSION=1 \n"); 
	fprintf(writefile, "TYPE=Scalar \n");
	
	for(int i=0; i<n_vertex; i++){
		fprintf( writefile, "%f\n", mWrite[i]);
	}
			
	fclose(writefile);
	return true;
}		


bool CCalculators::CAL_EMD(int n_BasicShape,int Test_N, std::vector<float> Basis_X, std::vector<float> SI_HIST, float* temp_EMD )
{
		
	for(int j=0;j<n_BasicShape;j++)
		{	
  		std::vector<nBin> IDEAL_SI_HIST(n_BasicShape, nBin(n_BasicShape,0));	                                 
			for(int i=0;i<n_BasicShape;i++)
			{	
				IDEAL_SI_HIST[i][i]=Test_N;	
			}
						
			TsSignature<string> srcSig(n_BasicShape, Basis_X, SI_HIST);
			TsSignature<string> snkSig(n_BasicShape, Basis_X, IDEAL_SI_HIST[j]);
		
			TsFlow flow[n_BasicShape];
			int flowVars = 0;
			double result = transportSimplex(&srcSig, &snkSig, DISTANCE, flow, &flowVars);
			
			//cout << "Total cost: " << result << endl;
			//cout << "Flows:" << endl;
		
			float temp1=0.0f;
			float temp2=0.0f;
			for (int k = 0; k < flowVars; k++)
			{ 
				//cout << basis_x[flow[i].from] << " to " << basis_x[flow[i].to] << " : " << flow[i].amount << endl;
				temp1 = temp1 + std::abs(Basis_X[flow[k].from]-Basis_X[flow[k].to])*flow[k].amount;
				temp2 = temp2 +flow[k].amount;	
			}
			
			temp_EMD[j] = temp1/temp2;
		}
	return true;
}


bool CCalculators::GoldenSection(int n_BasicShape,int Test_N, std::vector<float> Basis_X, std::vector<float> SI_HIST, float* temp_EMD, int STEP)
{ 
	typedef std::vector<float> IN;
	std::vector<IN> INIT_TEST_BIN(5,IN(n_BasicShape,0));
	
	std::vector<float> SLOPE(4,0);
	std::vector<nBin> IDEAL_SI_HIST(n_BasicShape, nBin(n_BasicShape,0));	                                 
	for(int i=0;i<n_BasicShape;i++)
	{	
		IDEAL_SI_HIST[i][i]=Test_N;	
	}
	
	// Initial Step
	for(int i=0; i<5; i++)
	{
		INIT_TEST_BIN[i] = IDEAL_SI_HIST[i*((n_BasicShape-1)/4)];
		
		TsSignature<string> srcSig(n_BasicShape, Basis_X, SI_HIST);
		TsSignature<string> snkSig(n_BasicShape, Basis_X, INIT_TEST_BIN[i]);
		
		TsFlow flow[n_BasicShape];
		int flowVars = 0;
		double result = transportSimplex(&srcSig, &snkSig, DISTANCE, flow, &flowVars);
		
		float temp1=0.0f;
		float temp2=0.0f;
		for (int k = 0; k < flowVars; k++)
		{ 
			temp1 = temp1 + std::abs(Basis_X[flow[k].from]-Basis_X[flow[k].to])*flow[k].amount;
			temp2 = temp2 +flow[k].amount;	
		}
		
		temp_EMD[i] = temp1/temp2;
	}
	
	for(int i=0;i<4; i++)
	{
		SLOPE[i] = temp_EMD[i+1] - temp_EMD[i];
	}
	int N_BETWEEN =100;
	for(int i=0;i<3; i++)
	{
		if(SLOPE[i] * SLOPE[i+1] <= 0)
			{
				N_BETWEEN = i;
			}
	}
	if(N_BETWEEN ==100)
	{
    if(SLOPE[0]<0 && SLOPE[1]<0 && SLOPE[2]<0 && SLOPE[3]<0)
    {
     	N_BETWEEN = 2;
    }
    if(SLOPE[0]>0 && SLOPE[1]>0 && SLOPE[2]>0 && SLOPE[3]>0)
    {
      N_BETWEEN = 0;
    }
  }
  
  if(N_BETWEEN ==100)
		{
  		cout << "N_BETWEEN = 100 !!!!!" << endl;	
  	}
  
	// Loop Step
  int tmp_START_POINT = 0;
  for(int j=0;j<(STEP-2);j++)
  {
  		int NUM_BIN = ((n_BasicShape-1)/pow(2,j+1));
      tmp_START_POINT = tmp_START_POINT + (NUM_BIN/2)*(N_BETWEEN);
     
      for(int i=0; i<5; i++)
      {
          INIT_TEST_BIN[i] = IDEAL_SI_HIST[tmp_START_POINT+(i)*(NUM_BIN/4)];
      }
      N_BETWEEN = OPT_EMD(temp_EMD,N_BETWEEN, Basis_X,INIT_TEST_BIN[1],INIT_TEST_BIN[3],SI_HIST,n_BasicShape);
  }
 
	return true;
}

int CCalculators::OPT_EMD(float* temp_EMD, int START,std::vector<float> BASE,std::vector<float> HIST_A2, std::vector<float> HIST_A4,std::vector<float> HIST_B, int n_BasicShape)
{
	std::vector<float> SLOPE(4,0);
	
	float* mEMD = new float[5];	
	mEMD[0] = temp_EMD[START];
	mEMD[2] = temp_EMD[START+1];
	mEMD[4] = temp_EMD[START+2];
	
  TsSignature<string> srcSig(n_BasicShape, BASE, HIST_B);
	TsSignature<string> snkSigA2(n_BasicShape, BASE, HIST_A2);
		
	TsFlow flow[n_BasicShape];
	int flowVars = 0;
	double result = transportSimplex(&srcSig, &snkSigA2, DISTANCE, flow, &flowVars);
		
	float temp1=0.0f;
	float temp2=0.0f;
	for (int k = 0; k < flowVars; k++)
	{ 
		temp1 = temp1 + std::abs(BASE[flow[k].from]-BASE[flow[k].to])*flow[k].amount;
		temp2 = temp2 +flow[k].amount;	
	}
		
	mEMD[1] = temp1/temp2;
	

	TsSignature<string> snkSigA4(n_BasicShape, BASE, HIST_A4);
	flowVars = 0;
	result = transportSimplex(&srcSig, &snkSigA4, DISTANCE, flow, &flowVars);
		
	temp1=0.0f;
	temp2=0.0f;
	for (int k = 0; k < flowVars; k++)
	{ 
		temp1 = temp1 + std::abs(BASE[flow[k].from]-BASE[flow[k].to])*flow[k].amount;
		temp2 = temp2 +flow[k].amount;	
	}
		
	mEMD[3] = temp1/temp2;	
	
	for(int i=0;i<4; i++)
	{
		SLOPE[i] = mEMD[i+1] - mEMD[i];
	}
	int N_BETWEEN =100;
	for(int i=0;i<3; i++)
	{
		if(SLOPE[i] * SLOPE[i+1] <= 0)
			{
				N_BETWEEN = i;
			}
	}
	
	if(N_BETWEEN ==100)
	{
    if(SLOPE[0]<0 && SLOPE[1]<0 && SLOPE[2]<0 && SLOPE[3]<0)
    {
     	N_BETWEEN = 2;
    }
    if(SLOPE[0]>0 && SLOPE[1]>0 && SLOPE[2]>0 && SLOPE[3]>0)
    {
      N_BETWEEN = 0;
    }
  }
  
  for(int i=0; i<5; i++)
  {
  	temp_EMD[i] = mEMD[i];
  }
  
  delete [] mEMD;
	mEMD=NULL;

	return N_BETWEEN;	
}		
	
double CCalculators::Cal_SA(Mesh* mesh, float* SA)
{
	double SUM_SA;
	const int n_vertex = mesh->nVertex();
	std::vector<float> n1;
	std::vector<float> n2;
	std::vector<float> pt;
		
	float edge1,edge2,edge3;
	
	for (int i = 0; i < n_vertex; i++)
	{
		pt.clear();
		const int *neighbor = mesh->vertex(i)->list();
		const int n_nbr = mesh->vertex(i)->nNeighbor();
		float temp =0.0f;
		
		for(int j=0; j<3; j++)
		{	
			pt.push_back(mesh->vertex(i)->fv()[j]);
			//cout << pt[j] << endl;
		}	
		
		for (int j = 0; j <  n_nbr-1; j++)
		{
			n1.clear();
			n2.clear();
			for(int k=0;k<3;k++)
			{
				n1.push_back(mesh->vertex(neighbor[j])->fv()[k]);
				n2.push_back(mesh->vertex(neighbor[j+1])->fv()[k]);
				//cout << n1[k] << endl;
			}
			
			edge1 = Distance_TwoPoint(pt,n1);
			edge2 = Distance_TwoPoint(pt,n2);
			edge3 = Distance_TwoPoint(n1,n2);
			temp = temp + Cal_TRI(edge1, edge2, edge3);
		}
		SA[i] = temp/n_nbr;
		SUM_SA += SA[i];
	}
	
	return SUM_SA;	
}

float CCalculators::Distance_TwoPoint(std::vector<float> pt1, std::vector<float> pt2)
{
	float Dist=0.0;
	Dist = sqrt( (pt1[0]-pt2[0])*(pt1[0]-pt2[0]) + (pt1[1]-pt2[1])*(pt1[1]-pt2[1]) + (pt1[2]-pt2[2])*(pt1[2]-pt2[2]) );
	
	return Dist;
}

float CCalculators::Cal_TRI(float edge1,float edge2,float edge3)
{
	float S = 0.0f;
	float m = 0.0f;
	m= (edge1+edge2+edge3)/2;
  S = sqrt(m*(m-edge1)+m*(m-edge2) + m*(m-edge3));
	return S;
}	         
        
