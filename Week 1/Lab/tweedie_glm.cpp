
// Space time
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y_i);

  PARAMETER( log_mu );
  PARAMETER( log_sd );
  PARAMETER( finv_power );
  PARAMETER_VECTOR( eps_j );
  
  Type jnll = 0; 
  Type sd = exp(log_sd);
  //Type power = 1.0 + invlogit(finv_power);
  Type power = 1.0 + exp(finv_power) / (1.0 + exp(finv_power));
  
  for( int i=0; i<Y_i.size(); i++){
     // dtweedie( Y_i(i), mean, sd, power, true );
     // where:
     // mean > 0
     // sd > 0
     // 1 < power < 2
     jnll -= dtweedie(Y_i(i), exp(log_mu), sd, power, true);
  }
  
  return jnll;
}
