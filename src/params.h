
namespace HSparams{
extern bool weakContribution ;
extern bool rescatter ;
extern bool shear ;
//extern double Temp, mu_b, mu_q, mu_s ;
extern int NEVENTS ;
extern double NBINS, QMAX ;
extern double dx, dy, deta ;
extern double ecrit ;

// ---- rooutines ----
void readParams(char* filename) ;
void printParameters() ;
}
