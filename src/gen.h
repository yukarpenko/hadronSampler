class TRandom3 ;
class DatabasePDG2;

namespace gen{
class Particle ;
//typedef std::vector<Particle*> ParticleList ; // TODO in far future
// data
extern DatabasePDG2 *database ;
extern TRandom3 *rnd ;
extern Particle ***pList ; // particle arrays
extern int *npart ;
const int NPartBuf = 100000; // dimension of particle buffer for each event

struct element {
 double tau, x, y, eta ;
 double u[4] ;
 double dsigma[4] ;
 double T, mub, muq, mus ;
 double pi[10] ;
 double Pi ;
} ;

// functions
void init();
void addElement(element elem) ;
int generate() ;
}

