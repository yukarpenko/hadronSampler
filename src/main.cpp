#include <omp.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TH1.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <TF1.h>
#include <sstream>
#include <TUUID.h>

#include "DatabasePDG2.h"
#include "gen.h"
#include "cascade.h"
#include "tree.h"
#include "params.h"

using namespace std ;
int getNlines(const char *filename) ;
void readCommandLine(int argc, char** argv) ;
int Nparticipants ;
string surfaceInput, eventOutput;

using HSparams::NEVENTS ;

int ranseed ;

extern "C"{
  void getranseedcpp_(int *seed)
  {
    *seed = ranseed ;
  }
}

// ########## MAIN block ##################

int main(int argc, char **argv)
{
 // command-line parameters
 readCommandLine(argc, argv) ;
 HSparams::printParameters() ;
 time_t time0 ;
 time(&time0) ;

 TRandom3* random3 = new TRandom3();
 random3->SetSeed(0);
 cout<<"hadronSampler: random seed = "<<random3->GetSeed()<<endl ;
 ranseed = random3->GetSeed();
 gen::rnd = random3 ;
  
//========= particle database init
	DatabasePDG2 *database = new DatabasePDG2("Tb/ptl3.data","Tb/dky3.mar.data");
	database->LoadData();
//	database->SetMassRange(0.01, 10.0); //-------without PHOTONS
//	database->SetWidthRange(0., 10.0);
	database->SortParticlesByMass() ;
	database->CorrectBranching() ;
	database->DumpData() ;
	cout << " pion index = " << database->GetPionIndex() << endl ;
  gen::database = database ;
  
// ========== generator init
 gen::load(surfaceInput.c_str(),getNlines(surfaceInput.c_str())) ;
 //cout << "dfMax = " << gen::calcDFMax() << endl ;

 // ========== trees & files
 time_t start, end ;
 time(&start);

//============= main task
 string eventOutputDir = eventOutput;
 eventOutputDir.erase(eventOutputDir.rfind("/"));
 char sbuffer [250];
 sprintf(sbuffer,"mkdir -p %s",eventOutputDir) ;
 system(sbuffer);
 TFile *outputFile = new TFile(eventOutput.c_str(), "RECREATE"); 
 outputFile->cd();
 MyTree *treeIni = new MyTree("treeini") ;
 MyTree *treeFin = new MyTree("treefin") ;
 
 gen::generate() ; // one call for NEVENTS

 for(int iev=0; iev<NEVENTS; iev++){
 treeIni->fill(iev, Nparticipants) ;
 gen::urqmd(iev) ;
 treeFin->fill(iev, Nparticipants) ;
 } // end events loop
 outputFile->Write() ;
 outputFile->Close() ;

 cout << "event generation done\n" ;
 time(&end); float diff2 = difftime(end, start);
 cout<<"Execution time = "<<diff2<< " [sec]" << endl;
 return 0;
}


void readCommandLine(int argc, char** argv)
{
  if(argc==1){
  cout << "no CL params - exiting.\n"; exit(1) ;
 }
 else{
  for(int iarg=1; iarg<argc-1; iarg++){
   if(strcmp(argv[iarg],"-Npart")==0) Nparticipants = atoi(argv[iarg+1]);
   if(strcmp(argv[iarg],"-params")==0) HSparams::readParams(argv[iarg+1]);
   if(strcmp(argv[iarg],"-surface")==0) surfaceInput = argv[iarg+1];
   if(strcmp(argv[iarg],"-output")==0) eventOutput = argv[iarg+1];
  }
  cout << "hadronSampler: command line parameters are:\n";
  cout << "Npart  " << Nparticipants << endl;
  cout << "hydro surface:  " << surfaceInput << endl;
  cout << "ROOT output:  " << eventOutput << endl;
 }
}


// auxiliary function to get the number of lines
int getNlines(const char *filename)
{
  ifstream fin(filename) ;
  if(!fin) {cout<<"getNlines: cannot open file: "<<filename<<endl; exit(1) ; }
  string line ;
  int nlines = 0 ;
  while(fin.good()){
    getline(fin,line) ; nlines++ ;
  } ;
  fin.close() ;
  return nlines-1 ;
}
