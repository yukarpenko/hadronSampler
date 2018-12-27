class TTree ;

class MyTree{
 TTree *tree ;
 Float_t *X, *Y, *Z, *T, *Px, *Py, *Pz, *E, *Temp, *mub, *muq, *mus ;
 Int_t *Id, *MId, *LastColl, *NColl, *Origin ;
 Short_t *Chrg, *Bar, *Strg ;
 Int_t nfill ;
public:
 MyTree(char *name) ;
 void fill(int iev) ;
} ;
