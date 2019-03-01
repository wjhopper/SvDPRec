#include <Rcpp.h>
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix diffusion_sim(int N, double a, double v, double t0,
                            double z, double sv, double st0,
                            double sz = 0, double s = 1,
                            NumericVector crit = {-.5, 5}) {
  
  NumericMatrix sim_data(N, 3);
  colnames(sim_data) = CharacterVector::create("RT", "speeded_resp","delayed_resp");
  z = z * a; // convert relative starting point to absolute
  double dt = .001; // timestep size
  
  NumericVector NDT = runif(N, t0 - st0/2, t0 + st0/2); // non-decision_times
  NumericVector evidence = rnorm(N, v, sv); // Sample SDT evidence strengths \ drift rates
  NumericVector drifts = evidence * dt; // Sample sampled evidence scale to instantaneous drift
  
  s = sqrt(pow(s, 2) * dt); // scale drift coefficient to instantaneous s.d.
  
  if (sz > 0) {
    NumericVector start_points = runif(N, z-.5*sz, z+.5*sz);
  } else {
    NumericVector z_vector = z;
    NumericVector start_points = rep_len(z_vector, N);
  }
  
  return(sim_data);
}
