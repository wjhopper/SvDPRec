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
  RVector<double> NDT;
  RVector<double> start_points;
  double a;
  double v;
  double sv;
  double s;
  RVector<double> crit;
  double dt;
  dqrng::normal_distribution rnorm;
  
  // constructor
  Diffusion(Rcpp::NumericMatrix sim_data,
            Rcpp::NumericVector NDT,
            Rcpp::NumericVector start_points,
            double a, double v, double sv, double s,
            Rcpp::NumericVector crit, double dt
            ) :
    sim_data(sim_data), NDT(NDT), start_points(start_points), a(a), v(v), sv(sv), s(s), crit(crit), dt(dt)
    {}
  
  void operator()(std::size_t begin, std::size_t end) {
    pcg64 rng(42, end);

    for(std::size_t i = begin; i < end; i++) {
      double evidence = v + rnorm(rng)*sv;  // Sample SDT evidence strength \ drift rate
      double drift = evidence*dt; // Sample sampled evidence scale to instantaneous drift
      
      int step = 0;
      double pos = start_points[i];
      
      while(pos < a && pos > 0 ) {
        pos += drift + rnorm(rng)*s;
        step += 1;
      }
      sim_data(i, 0) = NDT[i] + step*dt; // column 1 = rt
      
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
  
  // allocate the output vector
  NumericMatrix sim_data(N, 3);
  colnames(sim_data) = CharacterVector::create("RT", "speeded_resp","delayed_resp");
  z = z * a; // convert relative starting point to absolute
  
  NumericVector NDT = runif(N, t0, t0 + st0); // non-decision_times
  NumericVector start_points;
  if (sz > 0) {
    start_points = runif(N, z-.5*sz, z+.5*sz);
  } else {
    // NumericVector z_vector = NumericVector::create(z);
    start_points = rep_len(NumericVector::create(z), N);
  }
  Diffusion diffusion(sim_data, NDT, start_points, a, v, sv, s, crit, dt);
  RcppParallel::parallelFor(0, N, diffusion);
  return sim_data;
}