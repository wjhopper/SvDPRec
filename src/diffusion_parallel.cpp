#include <Rcpp.h>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(dqrng, BH, sitmo)]]
#include <pcg_random.hpp>
#include <dqrng_distribution.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

using namespace RcppParallel;

struct Diffusion : public Worker
{   
  RMatrix<double> sim_data;
  double a;
  double v;
  double t0;
  double z;
  double sv;
  double sz;
  double st0;
  double s;
  RVector<double> crit;
  double dt;

  // constructor
  Diffusion(Rcpp::NumericMatrix sim_data,
            double a, double v, double t0, double z, double sv, double sz, double st0, double s,
            Rcpp::NumericVector crit, double dt
            ) :
    sim_data(sim_data), a(a), v(v), t0(t0), z(z), sv(sv), sz(sz), st0(st0), s(s), crit(crit), dt(dt)
    {}
  
  void operator()(std::size_t begin, std::size_t end) {
    pcg64 rng(42, end);
    
    dqrng::normal_distribution evidence_dist(v, sv);
    dqrng::normal_distribution noise(0, s);
    dqrng::uniform_distribution NDT(t0, t0 + st0);
    dqrng::uniform_distribution start_points(z-.5*sz, z+.5*sz);
    
    for(std::size_t i = begin; i < end; i++) {
      double evidence = evidence_dist(rng);  // Sample SDT evidence strength \ drift rate
      double drift = evidence*dt; // Sample sampled evidence scale to instantaneous drift
      
      double pos;
      if (sz > 0) {
        pos = start_points(rng);
      } else{
        pos = z;
      }
      
      int step = 0;
      while(pos < a && pos > 0 ) {
        pos += drift + noise(rng);
        step += 1;
      }
      sim_data(i, 0) = NDT(rng) + step*dt; // column 1 = rt
      
      if (pos >= a) {
        sim_data(i, 1) = 1; // column 2 = speeded response
        if (evidence > crit[0]) {
          sim_data(i, 2) = 1; // column 3 = delayed response
        }
      } else {
        if (evidence > crit[1]) {
          sim_data(i, 2) = 1; // column 3 = delayed response
        }
      }
      
    }
  }
};

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector diffusion_parallel(int N, double a, double v, double t0,
                                       double z, double sv, double st0,
                                       double sz = 0, double s = 1,
                                       NumericVector crit = NumericVector::create(-.5, 5)) {
  
  double dt = .001; // timestep size
  s = sqrt(pow(s, 2.0l) * dt); // scale drift coefficient to instantaneous s.d.
  z = z * a; // convert relative starting point to absolute
  
  // allocate the output vector
  NumericMatrix sim_data(N, 3);
  colnames(sim_data) = CharacterVector::create("RT", "speeded_resp","delayed_resp");

  Diffusion diffusion(sim_data, a, v, t0, z, sz, sv, st0, s, crit, dt);
  RcppParallel::parallelFor(0, N, diffusion);
  return sim_data;
}