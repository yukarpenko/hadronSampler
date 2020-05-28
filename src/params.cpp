#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include "params.h"

using namespace std ;

namespace HSparams{

bool weakContribution ;
bool rescatter ;
bool shear ;
//double Temp, mu_b, mu_q, mu_s ;
int NEVENTS ;
double NBINS, QMAX ;
double dx, dy, deta ;
double ecrit ;

// ############ reading and processing the parameters

void readParams(const char* filename)
{
	char parName [255], parValue [255] ;
	ifstream fin(filename) ;
	if(!fin.is_open()) {
  cout << "Hadron Sampler: cannot open " << filename << endl;
  exit(1) ;
 }
 // setting default values
 NEVENTS = 10;
 rescatter = 1;
 weakContribution = true;
 shear = false;
 ecrit = 0.5;
 // reading the input file
	while(fin.good()){
	 string line ;
	 getline(fin, line) ;
	 istringstream sline (line) ;
	 sline >> parName >> parValue ;
	 if     (strcmp(parName,"number_of_events")==0) NEVENTS = atoi(parValue) ;
	 else if(strcmp(parName,"rescatter")==0) rescatter = atoi(parValue) ;
	 else if(strcmp(parName,"weakContribution")==0) weakContribution = atoi(parValue) ;
   else if(strcmp(parName,"shear")==0) shear = atoi(parValue) ;
   else if(strcmp(parName,"ecrit")==0) ecrit = atof(parValue) ;
	 else if(parName[0]=='!') cout << "CCC " << sline.str() << endl ;
	 else cout << "UUU " << sline.str() << endl ;
	}
 deta=0.05 ; dx=dy=0.0 ; // TODO!
 fin.close();
}

void printParameters()
{
  cout << "====== hadron sampler parameters ======\n" ;
  cout << "numberOfEvents = " << NEVENTS << endl ;
  cout << "isRescatter = " << rescatter << endl ;
  cout << "weakContribution = " << weakContribution << endl ;
  cout << "shear_visc_on = " << shear << endl ;
  cout << "e_critical = " << ecrit << endl ;
  cout << "======= end parameters =======\n" ;
}

}
