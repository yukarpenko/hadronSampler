#include <vector>

class TTree ;

namespace gen{
class Particle ;
}

class MyTree{
 TTree *tree ;
 Float_t *X, *Y, *Z, *T, *Px, *Py, *Pz, *E ;
 Int_t *Id, *MId, *LastColl, *NColl, *Origin, *Status ;
 Short_t *Chrg, *Bar, *Strg ;
 Int_t nfill ;
public:
 MyTree(char *name) ;
 void clear(void) { nfill = 0; }
 void addMediumHadrons(int iev);
 void addJetParticles(std::vector<gen::Particle>& ptls);
 void fillTree(void);
} ;
