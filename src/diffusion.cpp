#include <Rcpp.h>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(dqrng)]]
#include <dqrng.h>
// [[Rcpp::depends(BH, sitmo)]]
#include <pcg_random.hpp>
#include <dqrng_distribution.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

pcg_extras::seed_seq_from<std::random_device> seed_source;
pcg64 rng(seed_source);

// [[Rcpp::export]]
void set_diffusion_SDT_seed(uint64_t new_seed) {
  rng.seed(new_seed);
}

// The following diffusion_SDT_sim function is intentially not exported to R.
// It exists in case the need for some backwards compatibility arises, or if in the future
// the functionallity to switch between different parallel and serial implementations arises.
// Currently, the recommended approach to switch between parallel and serial execution is to restrict
// the RcppParallel library to a single CPU core by calling RcppParallel::setThreadOptions(numThreads = 1)
// from within the R session when you need serial, and calling RcppParallel::setThreadOptions(numThreads = 'auto')
// when you want parallelism.

using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix diffusion_SDT_sim2(int N, double a, double v, double t0,
                                double z, double sz = 0, double sv = 0,
                                double st0 = 0, double s = 1,
                                NumericVector crit = NumericVector::create(0, 0)) {

  dqrng::dqRNGkind("pcg64");
  dqrng::dqset_seed(abs(rng()));
  rng.backstep(1);

  NumericMatrix sim_data(N, 3);
  colnames(sim_data) = CharacterVector::create("RT", "speeded_resp","delayed_resp");
  double dt = .001; // timestep size
  // 
  NumericVector evidence = dqrng::dqrnorm(N, v, sv); // Sample SDT evidence strengths \ drift rates
  NumericVector drifts = evidence * dt; // Sample sampled evidence scale to instantaneous drift
  // 
  s = sqrt(pow(s, 2.0l) * dt); // scale drift coefficient to instantaneous s.d.
  // 
  NumericVector start_points(no_init(N));
  if (sz > 0) {
    start_points = dqrng::dqrunif(N, z-.5*sz, z+.5*sz);
  } else {
    start_points = NumericVector(N, z);
  }
  
  NumericVector NDT(no_init(N)); // non-decision_times
  if (st0 > 0) {
    NDT = dqrng::dqrunif(N, t0, t0 + st0);
  } else {
    NDT = NumericVector(N, t0);
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

int nThreads() {
  // determine max number of threads
  std::size_t threads = tthread::thread::hardware_concurrency();
  char* numThreads = ::getenv("RCPP_PARALLEL_NUM_THREADS");
  if (numThreads != NULL) {
    int parsedThreads = ::atoi(numThreads);
    if (parsedThreads > 0)
      threads = parsedThreads;
  }
  return threads;
}
// To-do Notes: This following parallel implementation based on the RcppParallel framework
// currently uses a thread-local instantiation of the PCG pseudo-random number generator
// with a fixed seed. Currently, the only way to set the prng seed programmatically 
// would be adding a "seed" argument to the function interfacethat is exported to R, 
// and then passing this seed on to the constructor for the RcppParallel::Worker
// derived object. I do not like this idea, it clutters the API for the diffusion model and
// doesn't divide separate functionallity (seeding the rng) into separate functions. So, I am
// leaving it with a fixed seed for now.
//
// If the need to programmatically seed the prng arises in the future (e.g., from an R session),
// there are two possible modifications I can envision. First, use a global pcg64 instance,
// write a function which uses the .seed() method to set the seed, then use this global
// object in each thread, switching to a thread-specific prng stream with the set_stream() method.
// Alternatively, the seed could be a global variable, have a function to modify this global value,
// and then the thread-local pcg64 instances could use this global value as their seed.


using namespace RcppParallel;

struct Diffusion : public Worker
{   
  RMatrix<double> sim_data;
  double a;
  double v;
  double t0;
  double z;
  double sz;
  double sv;
  double st0;
  double s;
  RVector<double> crit;
  double dt;

  // constructor
  Diffusion(Rcpp::NumericMatrix sim_data,
            double a, double v, double t0, double z, double sz, double sv, double st0, double s,
            Rcpp::NumericVector crit, double dt
  ) :
    sim_data(sim_data), a(a), v(v), t0(t0), z(z), sz(sz), sv(sv), st0(st0), s(s), crit(crit), dt(dt)
  {}
  
  void operator()(std::size_t begin, std::size_t end) {
    rng.set_stream(end);

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

      if (st0 > 0) {
        sim_data(i, 0) = NDT(rng) + step*dt; // column 1 = rt
      } else {
        sim_data(i, 0) = t0 + step*dt; // column 1 = rt
      }

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
Rcpp::NumericVector diffusion_SDT(int N, double a, double v, double t0,
                                  double z, double sz = 0, double sv = 0, 
                                  double st0 = 0, double s = 1,
                                  NumericVector crit = NumericVector::create(0, 0)) {

  double dt = .001; // timestep size
  s = sqrt(pow(s, 2.0l) * dt); // scale drift coefficient to instantaneous s.d.

  // allocate the output vector
  NumericMatrix sim_data(N, 3);
  colnames(sim_data) = CharacterVector::create("RT", "speeded_resp","delayed_resp");
  
  Diffusion diffusion(sim_data, a, v, t0, z, sz, sv, st0, s, crit, dt);
  RcppParallel::parallelFor(0, N, diffusion, N/nThreads());
  
  sim_data.attr("class") = CharacterVector::create("diffusion_SDT","matrix");
  return sim_data;
}
