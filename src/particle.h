class ParticlePDG2 ;

namespace gen {
class Particle {
public:
  double px, py, pz, e ;
  double x, y, z, t ;
  ParticlePDG2 *def ;
  int mid, ncoll, lastcoll, origin ;
  int status;
  Particle(double X, double Y, double Z, double T,
    double Px, double Py, double Pz, double E, ParticlePDG2* def1, int Mid):
    px(Px), py(Py), pz(Pz), e(E), x(X), y(Y), z(Z), t(T), def(def1), mid(Mid),
    ncoll(0), lastcoll(0), origin(0), status(0) {} ;
  Particle(double X, double Y, double Z, double T,   // constructor with status argument
    double Px, double Py, double Pz, double E, ParticlePDG2* def1, int Mid, int Status):
    px(Px), py(Py), pz(Pz), e(E), x(X), y(Y), z(Z), t(T), def(def1), mid(Mid),
    ncoll(0), lastcoll(0), origin(0), status(Status) {} ;
  ~Particle() {} ;
} ;
}
