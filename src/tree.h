class TTree ;

class MyTree{
 TTree *tree ;
 Float_t *X, *Y, *Z, *T, *Px, *Py, *Pz, *E ;
 Int_t *Id, *MId, *LastColl, *NColl, *Origin ;
 Short_t *Chrg, *Bar, *Strg ;
 Int_t nfill ;
 Int_t np ;
public:
 MyTree(char *name) ;
 void fill(int iev, int Nparticipants) ;
} ;
