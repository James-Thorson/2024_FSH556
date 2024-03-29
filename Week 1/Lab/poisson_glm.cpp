
// Space time
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y_i);
  PARAMETER( log_mu );
  Type jnll = 0; 
  for( int i=0; i<Y_i.size(); i++){
     jnll -= dpois(Y_i(i), exp(log_mu), true);
  }
  return jnll;
}
