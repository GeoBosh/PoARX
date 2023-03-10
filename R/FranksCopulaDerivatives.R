# first derivative of generator w.r.t. theta
.d1dtheta_iPsiFrank <- function(t, theta){
  #### 1/(exp(theta)-1) - t/(exp(theta*t)-1)
  1/expm1(theta) - t/expm1(theta*t)
}

# first derivative of generator w.r.t. t
.d1dt_iPsiFrank <- function(t, theta){
  #### -theta/(exp(theta*t)-1)
  - theta/expm1(theta*t)
}

# second derivative of generator w.r.t. theta (twice)
.d2dtheta2_iPsiFrank <- function(t, theta){
  #### t^2*exp(theta*t)/(exp(theta*t)-1)^2 - exp(theta)/(exp(theta)-1)^2
  t^2*(exp(theta*t)) / (expm1(theta*t)^2) -
    (exp(theta)) / ( expm1(theta)^2)
}

# second derivative of generator w.r.t. theta and t
.d2dthetat_iPsiFrank <- function(t, theta){
  #### (1 - exp(theta*t) + theta*t*exp(theta*t))/(exp(theta*t)-1)^2
  etht <- expm1(theta*t)
  (theta*t*(exp(theta*t)) - etht) / (etht^2)
}

# first derivative of inverse generator w.r.t. theta
.d1dtheta_psiFrank <- function(t, theta){
  #### 1/theta^2 * log(1 + exp(-t)(exp(-theta)-1)) +
  ####   exp(- theta - t) / (theta*(1 + exp(-t)(exp(-theta)-1)))
  bit <- expm1(-theta)*exp(-t)
  1/theta^2* log1p(bit) + (exp(-t-theta))/(theta*(1+bit))
}

# first derivative of inverse generator w.r.t. t
.d1dt_psiFrank <- function(t, theta){
  #### exp(-t)(exp(-theta)-1)/(theta*(1 + exp(-t)(exp(-theta)-1))
  bit <- expm1(-theta)*exp(-t)
  bit / (theta * (1 + bit))
}

# second derivative of inverse generator w.r.t. theta (twice)
.d2dtheta2_psiFrank <- function(t, theta){
  #### exp(- theta - t)*(exp(-t)-1) / (theta*(1 + exp(-t)(exp(-theta)-1))^2)-
  ####   2exp(- theta - t)/(theta^2*(1 + exp(-t)(exp(-theta)-1))) -
  ####   2/theta^3*log(1 + exp(-t)(exp(-theta)-1))
  bit <- expm1(-theta)*exp(-t)
  (exp(-t-theta))*expm1(-t) / (theta*(1+bit)^2) -
    2 * (exp(-t-theta)) / (theta^2*(1+bit)) -
    2 * log1p(bit) / theta^3
}

# second derivative of inverse generator w.r.t. theta and t
.d2dthetat_psiFrank <- function(t, theta){
  #### exp(- theta - t)*(exp(-t)-1) / (theta*(1 + exp(-t)(exp(-theta)-1))^2)-
  ####   2exp(- theta - t)/(theta^2*(1 + exp(-t)(exp(-theta)-1))) -
  ####   2/theta^3*log(1 + exp(-t)(exp(-theta)-1))
  bit <- 1+expm1(-theta)*exp(-t)

  - expm1(-theta)*(exp(-t)) / (theta^2*bit) -
    (exp(-t-theta)) / (theta*bit^2)
}

# second derivative of inverse generator w.r.t. t (twice)
.d2dt2_psiFrank <- function(t, theta){
  #### ((exp(-t)(exp(-theta)-1))^2 -
  ####   (1 + exp(-t)(exp(-theta)-1))*(exp(-t)(exp(-theta)-1))) /
  ####     (theta*(1 + exp(-t)(exp(-theta)-1)))^2
  bit <- expm1(-theta)*exp(-t)
  - bit / (theta * (1+bit)^2)
}
