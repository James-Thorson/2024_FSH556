
// Space time
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y_i);
  DATA_IVECTOR(j_i);

  PARAMETER( log_mu );
  PARAMETER( log_sd_indiv );
  PARAMETER_VECTOR( eps_j );
  
  Type jnll = 0; 
  
  // Distribution for random effect
  for( int j=0; j<eps_j.size(); j++ ){
     jnll -= dnorm( eps_j(j), Type(0.0), exp(log_sd_indiv), true);
  }
  
  // Conditional probability of data
  vector<Type> lambda_i( Y_i.size() );
  vector<Type> eps_i( Y_i.size() );
  for( int i=0; i<Y_i.size(); i++){
     eps_i(i) = eps_j(j_i(i));
     lambda_i(i) = exp( log_mu + eps_i(i) );
     jnll -= dpois(Y_i(i), lambda_i(i), true);
  }
  
  return jnll;
}
