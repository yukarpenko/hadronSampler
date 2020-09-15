class TRandom3 ;
class DatabasePDG2;
class Particle ;

namespace gen{
//typedef std::vector<Particle*> ParticleList ; // TODO in far future
// data
extern DatabasePDG2 *database ;
extern TRandom3 *rnd ;
extern Particle ***pList ; // particle arrays
extern int *npart ;
const int NPartBuf = 150000; // dimension of particle buffer for each event

// functions
void load(const char *filename, int N) ;
int generate() ;
}

