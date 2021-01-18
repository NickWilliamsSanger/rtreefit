functions{

    real x_to_t2(real[] x, int[] parentidx, int[] xidx,real[] tip_min_age, real[] t, int i);
    //real[] xx_to_t2(real[] x, int[] parentidx, int[] xidx,real[] tip_min_age, real[] t, int i);
    //TODO (See above). The following is a transcription of the C code that supported pass-by-reference
    //pass-by-ref is not supported in stan - so we can set each t multiple times. Should be possible to fix
    //by passing back a vector.
    real x_to_t2(real[] x, int[] parentidx, int[] xidx,real[] tip_min_age, real[] t, int i){
      real factor;
      if(t[i]>=0){
        return t[i];
      }
      if(parentidx[i]<1){
        if(xidx[i]<1){
          //tip is directly ancestral to root..
          return tip_min_age[i];
        }else{
          //t[i]=x[xidx[i]]*tip_min_age[i];
          return x[xidx[i]]*tip_min_age[i];
        }
      }else{
        real tmp=0;
        int k=i;
        while(parentidx[k]>=1){
          k=parentidx[k];
          tmp+=x_to_t2(x,parentidx,xidx,tip_min_age,t,k);
        }
        if(xidx[i]>=1){
          //interior
          factor=tip_min_age[i]-tmp;
          //t[i]=x[xidx[i]]*factor;
          return x[xidx[i]]*factor;
        }else{
          //t[i]=tip_min_age[i]-tmp;
          return tip_min_age[i]-tmp;
        }
        //return t[i];
      }
    }
    //Need to make sure all indices are 1+
    //Transforms stick-breaking fractions, x, to branch durations, t
      vector x_to_t(real[] x, int[] parentidx, int[] xidx, real[] tip_min_age, int n){
        int N=n;
        real t[N];
        for(i in 1:N){
          t[i]=-1;
        }
        for(i  in 1:N){
          t[i]=x_to_t2(x,parentidx,xidx,tip_min_age,t,i);
        }
        return to_vector(t);
      }

      vector logisticMean(real L,real k,real midpoint,vector a,vector b){
        return ((L/k)*(log(1+exp(k*(b-midpoint)))-log(1+exp(k*(a-midpoint))))) ./ (b-a);
      }

      vector getExtraLambdaRates(vector t,int[] parentidx,int N){
        real extralambda=149.1968140;
        real kg=-50.0161248;
        real mp=0.2246393;
        real ct=0;//current time
        int k;
        vector[N] t0;

        for(i in 1:N){
          t0[i]=-1;
        }
        for(i in 1:N){
          k=parentidx[i];
          ct=0.0;
          while(k>0){
            if(t0[k]>0){
              ct=ct+t0[k]+t[k];
              k=-1;
            }else{
              ct=ct+t[k];
              k=parentidx[k];
            }

          }
          t0[i]=ct;
        }
        for(i in 1:N){
          if (t0[i]<0){
          reject("x must not be negative; found x=", t0[i]);
          }

        }
        return logisticMean(extralambda,kg,mp,t0,t0+t);
      }



  }

data{
  int N; //num branches
  int NINT; //num internal branches
  int NLAMBDA;//number of rates to infer
  int parentidx[N];//Index of parent branch (-1 for root)
  int xidx[N];
  real tip_min_age[N];
  int rates[N];
  int ratesp[N];
  vector[N] s; //sensitivity
  int m[N];
  int nh[NINT];
  vector[NINT] q;
  vector[NINT] concentration;
  int  idxcrossover[NLAMBDA-1];
  real lambda_est;
  real early_growth_model_on;
}

parameters {
  real<lower=0.0001,upper=0.9999> x[NINT];
  //vector<lower=0.001,upper=0.999>[N] S;
  vector<lower=1,upper=200>[NLAMBDA] lambda;
  real<lower=0,upper=1> x0[NLAMBDA-1]; //fractional crossover..
  real<lower=0> k;  //nb 1/overdispersion
}

model {
  vector[N] t;
  vector[N] t0;
  vector[N] lambda_per_branch;
  vector[N] tmp;
  vector[N] x0v;
  x0v=rep_vector(0.0,N);

  //real t0=0.0;
  x ~ beta(concentration .* q ./ (1-q),concentration);
  //S ~ beta(100,100*(1-s) ./ s);
  x0 ~ uniform(0,1);
  k ~ normal(0,1);
  for(i in 1:(NLAMBDA-1)){
    x0v[idxcrossover[i]]=x0[i];
  }

  lambda ~ normal(lambda_est,0.25*lambda_est);
  t=x_to_t(x, parentidx, xidx, tip_min_age, N);
  lambda_per_branch=(1-x0v) .* lambda[rates]+x0v .* lambda[ratesp]+early_growth_model_on*getExtraLambdaRates(t,parentidx,N);
  //under thinning the negative binomial is again negative binomial..
  //See ""Puig, Pedro, and Jordi Valero. "Characterization of count data distributions involving additivity and binomial subsampling." Bernoulli (2007): 544-555.""
  // key point is that mean -> s*mean and that (variance-mean)/mean**2 is invariant under subsampling
  // However under changes in time scale,  (variance-mean)/mean**2 varies linearly with time so we need to
  // scale phi by the underlying mean (t*lambda)..
  m ~ neg_binomial_2(t .* lambda_per_branch .* s , t .* lambda_per_branch/k);
}

generated quantities {
  vector[N] ta;
  ta=x_to_t(x, parentidx, xidx, tip_min_age, N);
}
