#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat center_rows(const arma::mat & G){
  int N;
  arma::mat G1;
  arma::mat G2;
  
  G1 = G.t();
  
  N = G1.n_rows;
  G2 = G1 - arma::repmat(arma::mean(G1), N, 1);
  
  return G2.t();
}

arma::mat scale_rows(const arma::mat & G){
  int N;
  arma::mat G1;
  arma::mat G2;
  
  G1 = G.t();
  
  N = G1.n_rows;
  G2 = G1 / arma::repmat(arma::stddev(G1), N, 1);
  
  return G2.t();
}

double mentapprc(const arma::rowvec & x){
  // define variables
  int n;
  arma::rowvec x_;
  double xstd, k1, k2, gamma, gaussianEntropy, negentropy, entropy;
  
  // compute size of vector
  n = x.size();
  
  // standardize vector
  x_ = scale_rows(center_rows(x));
  
  // constants we need
  xstd = arma::stddev(x);
  k1 = 36 / (8 * sqrt(3) - 9);
  gamma = 0.37457;
  k2 = 79.047;
  gaussianEntropy = log(2 * M_PI) / 2 + 0.5;
  
  // this is negentropy
  negentropy = k2 * pow((arma::mean(arma::log(arma::cosh(x_))) - gamma), 2)
    + k1 * pow(arma::mean(x_ % arma::exp(- pow(x_, 2) / 2)), 2);
  
  // this is entropy
  entropy = gaussianEntropy - negentropy + log(xstd);
  
  // return result
  return entropy;
}

double pwlingc(const arma::mat & X) {
  // define variables
  int p, n;
  arma::mat X_;
  arma::mat C;
  arma::rowvec res1;
  arma::rowvec res2;
  double LR;
  
  // compute size of matrix
  p = X.n_rows;
  n = X.n_cols;
  
  // standardize matrix
  X_ = scale_rows(center_rows(X));
  
  // compute covariance matrix
  C = cov(X_.t());
  
  
  // general entropy-based method, for variables of any distribution
  res1 = X_.row(0) - C(0, 1) * X_.row(1);
  res2 = X_.row(1) - C(1, 0) * X_.row(0);
  LR = mentapprc(X_.row(0)) - mentapprc(X_.row(1))
    - mentapprc(res1) + mentapprc(res2);
  
  // Rcpp::Rcout << res1 << "\n"
  // << res2 << "\n";
  
  // return result
  return pow(std::min(0.0, LR), 2);
}

arma::vec mysetdiff(const arma::vec & x, const arma::vec & y) {
  std::vector<double> diff;
  
  std::set_difference(x.begin(), x.end(), y.begin(), y.end(),
                      std::inserter(diff, diff.begin()));
  
  return arma::vec(diff);
  
}

arma::mat bind_row_vectors(const arma::rowvec & A, const arma::rowvec & B){
  arma::mat out = join_vert(A, B); // same as join_cols
  return out;
}

int findindexc(const arma::mat & X,
               const arma::vec & candidates,
               const arma::vec & U_K){
  // define variables
  int i, j, n, p;
  double minT, J;
  
  n = X.n_cols;
  p = X.n_rows;
  arma::vec T_vec(p);
  arma::mat X_temp(2, n);
  
  // initialize vars
  T_vec.fill(NA_REAL);
  
  // Rcpp::Rcout << T_vec;
  
  minT = -1;
  
  // loop through vars
  for (j = 0; j < candidates.size(); j++){
    
    int candidate_j = candidates[j];
    
    if (minT == -1){
      
      T_vec[candidate_j] = 0;
      
      for (i = 0; i < U_K.size(); i++){
        
        int U_K_i = U_K[i];
        
        if (U_K_i == candidate_j){
          
          continue;
          
        } else {
          
          X_temp = bind_row_vectors(X.row(U_K_i), X.row(candidate_j));
          
          // Rcpp::Rcout << "Number of cols: " << X_temp.n_cols
          //             << ", number of rows: " << X_temp.n_rows << "\n";
          
          J = pwlingc(X_temp);
          T_vec[candidate_j] += J;
          
          // Rcpp::Rcout << "Value for J: " << J
          //             << ", value for T_vec: " << T_vec[j] << "\n";
        }
      }
      
      minT = T_vec[candidate_j];
      
    } else {
      
      T_vec[candidate_j] = 0;
      
      for (i = 0; i < U_K.size(); i++){
        
        int U_K_i = U_K[i];
        
        if (U_K_i == candidate_j){
          
          continue;
          
        } else {
          
          X_temp = bind_row_vectors(X.row(U_K_i), X.row(candidate_j));
          
          // Rcpp::Rcout << "Number of cols: " << X_temp.n_cols
          //             << ", number of rows: " << X_temp.n_rows << "\n";
          
          J = pwlingc(X_temp);
          T_vec[candidate_j] += J;
          
          // Rcpp::Rcout << "Value for J: " << J
          //             << ", value for T_vec: " << T_vec[j] << "\n";
          
          if (T_vec[candidate_j] > minT){
            
            T_vec[candidate_j] = std::numeric_limits<double>::infinity();
            break;
            
          }
        }
      }
      
      minT = std::min(T_vec[candidate_j], minT);
      
    }
  }
  
  // Rcpp::Rcout << T_vec;
  
  return index_min(T_vec) + 1;
}

arma::cube computeRc(const arma::mat & X,
                     const arma::vec & candidates,
                     const arma::vec & U_K,
                     const arma::mat & M){
  // define variables
  int i, j, n, p;
  double minT, J;
  
  n = X.n_cols;
  p = X.n_rows;
  arma::cube R = arma::cube(p, n, p, arma::fill::zeros);
  
  // compute cov matrix
  arma::mat Cov = arma::cov(X.t());
  
  for (j = 0; j < candidates.size(); j++){
    
    int candidate_j = candidates[j];
    
    for (i = 0; i < U_K.size(); i++){
      
      int U_K_i = U_K[i];
      
      if (U_K_i == candidate_j){
        
        continue;
        
      } else if (M.at(U_K_i, candidate_j) == 0){
        // Rcpp::Rcout << i << j << "\n";
        R(arma::span(U_K_i), arma::span::all, arma::span(candidate_j)) = X.row(U_K_i);
      } else {
        R(arma::span(U_K_i), arma::span::all, arma::span(candidate_j)) =
          X.row(U_K_i) - Cov.at(U_K_i, candidate_j) / Cov.at(candidate_j, candidate_j)
        * X.row(candidate_j);
      }
    }
  }
  
  return R;
}


// [[Rcpp::export]]
std::vector<double> direct_lingam_c(const arma::mat & X){
  // define variables
  int n, p, index;
  
  n = X.n_rows;
  p = X.n_cols;
  
  arma::mat X_(p, n);
  arma::mat M(p, p);
  arma::rowvec K(p);
  arma::vec candidates;
  arma::cube R(p, n, p);
  
  // transpose matrix
  X_ = X.t();
  
  // center variables
  X_ = center_rows(X_);
  
  // prepare matrix M
  M.fill(-1);
  M.diag().fill(0);
  
  // prepare vector K and U_K
  K = arma::zeros<arma::rowvec>(p);
  candidates = arma::linspace(0, p - 1, p);
  
  // loop through variables
  for (int m = 0; m < p; m++){
    
    // compute R, i.e., matrix containing residuals
    R = computeRc(X_, candidates, candidates, M);
    
    // find exogenous variables
    if (candidates.size() == 1){
      index = candidates[0];
    } else {
      index = findindexc(X_, candidates, candidates) - 1;
    }
    
    // update causal order
    K[m] = index;
    
    
    // update objects
    M.col(index).fill(NA_REAL);
    M.row(index).fill(NA_REAL);
    X_ = R(arma::span::all, arma::span::all, arma::span(index));
    
    // update candidates
    // Rcpp::Rcout << "m = " << m << "\n";
    // Rcpp::Rcout << "one-by-one = " << arma::vec(1).fill(index) << "\n";
    // Rcpp::Rcout << "index = " << index << "\n";
    // Rcpp::Rcout << "candidates = " <<
    // mysetdiff(candidates, arma::vec(1).fill(index)) << "\n";
    candidates = mysetdiff(candidates, arma::vec(1).fill(index));
    
  }
  
  // return causal order
  K = K + 1;
  
  return arma::conv_to<std::vector<double>>::from(K);
  
}