#include <Rcpp.h>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(dqrng)]]
#include <dqrng.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix diffusion_SDT_sim(int N, double a, double v, double t0,
                                double z, double sv, double st0,
                                double sz = 0, double s = 1,
                                NumericVector crit = NumericVector::create(-.5, 5)) {

  dqrng::dqRNGkind("pcg64");
  dqrng::dqset_seed(42);

  NumericMatrix sim_data(N, 3);
  colnames(sim_data) = CharacterVector::create("RT", "speeded_resp","delayed_resp");
  z = z * a; // convert relative starting point to absolute
  double dt = .001; // timestep size
  // 
  NumericVector NDT = dqrng::dqrunif(N, t0, t0 + st0); // non-decision_times
  NumericVector evidence = dqrng::dqrnorm(N, v, sv); // Sample SDT evidence strengths \ drift rates
  NumericVector drifts = evidence * dt; // Sample sampled evidence scale to instantaneous drift
  // 
  s = sqrt(pow(s, 2.0l) * dt); // scale drift coefficient to instantaneous s.d.
  // 
  NumericVector start_points;
  if (sz > 0) {
    start_points = dqrng::dqrunif(N, z-.5*sz, z+.5*sz);
  } else {
    start_points = NumericVector(N, z);
  }
  
    for(int i = 0; i < N; i++) {
      int step = 0;
      double pos = start_points[i];
      double v_inst = drifts[i];
    
      while(!(pos > a || pos < 0 )) {
        step += 1;
        pos +=  dqrng::rnorm(v_inst, s);
      }
    
      if (pos >= a){
        sim_data(i, 1) = 1; // column 2 = speeded response
      }
  
      sim_data(i, 0) = step*dt; // column 1 = rt
    }
    
    sim_data(_, 0) =  NDT + sim_data(_, 0);
  
    LogicalVector above_crit_o = evidence > crit[0];
    LogicalVector above_crit_n = evidence > crit[1];
    LogicalVector said_old = sim_data(_, 1) == 1;
  
    for (int j = 0; j < sim_data.nrow(); j++) {
      if ((above_crit_o[j] && said_old[j]) || (above_crit_n[j] && !said_old[j])) {
        sim_data(j, 2) = 1;
      }
    }

  return(sim_data);
}
