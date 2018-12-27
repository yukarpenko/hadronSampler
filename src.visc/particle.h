class ParticlePDG2 ;

class Particle {
public:
  double px, py, pz, e ;
  double x, y, z, t ;
  double Temp, mub, muq, mus; // to keep track of T/mu_i where the particle was produced
  ParticlePDG2 *def ;
  int mid, ncoll, lastcoll, origin ;
  Particle(double X, double Y, double Z, double T,
    double Px, double Py, double Pz, double E,
    double Temp, double mub, double muq, double mus, ParticlePDG2* def1, int Mid):
    px(Px), py(Py), pz(Pz), e(E), x(X), y(Y), z(Z), t(T),
    Temp(Temp), mub(mub), muq(muq), mus(mus), def(def1), mid(Mid),
    ncoll(0), lastcoll(0), origin(0) {} ;
  ~Particle() {} ;
} ;
