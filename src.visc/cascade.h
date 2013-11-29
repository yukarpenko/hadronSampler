
extern"C" {
 int geteposcode_(int *code);
 void pdg2id_(int *ityp, int *iso3, int *pdgid) ;
 void urqmdmain_();

 void cxxinit_(int* index, int* id, float* x, float* y, float* z, float* t, float* px, float* py, float* pz, float* E, float* mass);
 void cxxninit_(int *np) ;
 void cxxnfinal_(int *np) ;
 void cxxfinal_(int* index, int* id, float* x, float* y, float* z, float* t, float* px, float* py, float* pz, float* E, float* mass);
}

void decay(int pid, double mom[4], int& nprod, int ppid [], double pmom [][4]) ;

namespace gen{
void urqmd(int ievent) ;
}
