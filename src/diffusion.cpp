#include <Rcpp.h>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppZiggurat)]]
#include <Ziggurat.h>

using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::export]]
void zigg_seed(unsigned long int s) {
  zigg.setSeed(s);
  return;
}

// [[Rcpp::export]]
NumericMatrix diffusion_SDT_sim(int N, double a, double v, double t0,
                                double z, double sv, double st0,
                                double sz = 0, double s = 1,
                                NumericVector crit = NumericVector::create(-.5, 5)) {
  
  NumericMatrix sim_data(N, 3);
  colnames(sim_data) = CharacterVector::create("RT", "speeded_resp","delayed_resp");
  z = z * a; // convert relative starting point to absolute
  double dt = .001; // timestep size
  // 
  NumericVector NDT = runif(N, t0, t0 + st0); // non-decision_times
  NumericVector evidence = rnorm(N, v, sv); // Sample SDT evidence strengths \ drift rates
  NumericVector drifts = evidence * dt; // Sample sampled evidence scale to instantaneous drift
  // 
  s = sqrt(pow(s, 2) * dt); // scale drift coefficient to instantaneous s.d.
  // 
  NumericVector start_points;
  if (sz > 0) {
    start_points = runif(N, z-.5*sz, z+.5*sz);
  } else {
    NumericVector z_vector = NumericVector::create(z);
    start_points = rep_len(z_vector, N);
  }
  
    for(int i = 0; i < N; i++) {
      int step = 0;
      double pos = start_points[i];
      double v_inst = drifts[i];
    
      while(!(pos > a || pos < 0 )) {
        step += 1;
        pos +=  v_inst + zigg.norm()*s;
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
