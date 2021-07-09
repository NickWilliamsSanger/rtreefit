require("phytools")
require("rstan")
#' Fit tree to infer branch timings and mutation rate based on observed mutation count and sensitivity
#'
#' @export
#' @param tree Augmented ape phylo object. A list that should contain agedf - a dataframe with tip.label and age specifying the terminal timepoint for each tip.
#' @param switch_nodes Nodes where a change rate occurs
#' @param xcross  Currently ignored.  Method assumes uniform prior for location of switch event on branch.  Plan to reintroduce this to fix distance.
#' @param b_pool_rates Boolean. Whether to pool the mutant clades specified by switch_nodes and estimate as single mutant mutation rate
#' @param niter  Number of iterations per chain in stan inference
#' @param model nb_tree or poisson_tree
#' @param early_growth_model_on  Numeric. Scaling for early growth model.  If 0 then the model is switched off.
#' @param stan_control  List. Control list to be passed into rstan::sampling.
#' @return A list
#'
fit_tree=function(tree,
                  switch_nodes,
                  xcross=NA,
                  b_pool_rates=FALSE,
                  niter=20000,
                  cores=3,
                  split_nodes=switch_nodes,
                  model="nb_tree",
                  early_growth_model_on=1.0,
                  stan_control=list(adapt_delta=0.95)
                  ){
  if(is.null(tree$sensitivity)){
    warning("No sensitivity supplied: assuming 99%")
    tree$sensitivity=rep(0.99,length(tree$edge.length))
  }
  if(is.null(tree$agedf)){
    stop("Please supply dataframe: agedf")
  }
  pt=list(tree=tree,agedf=tree$agedf,switch_node=switch_nodes,xcross=xcross)
  fitres=nbfit_tree_stan_mcmc(tree = pt$tree,
                              agedf = pt$agedf,
                              rate_switch_nodes = pt$switch_node,
                              xcross=pt$xcross,
                              niter = niter,
                              b_pool_rates = b_pool_rates,
                              cores=cores,
                              split_nodes = split_nodes,
                              model = model,
                              early_growth_model_on=early_growth_model_on,
                              stan_control = stan_control)

  fitres$agedf=pt$agedf
  fitres
}

nbfit_tree_stan_mcmc=function(tree,
         agedf,
         rate_switch_nodes,
         xcross=rep(0.99,length(rate_switch_nodes)),
         b_pool_rates=TRUE,
         niter=50000,
         cores=3,
         b_debug=FALSE,model="nb_tree",
         odf=-1,
         split_nodes=rate_switch_nodes,
         concentration=1,
         early_growth_model_on=1.0,
         stan_control=list(adapt_delta=0.95)
){
  if(b_debug){
    browser()
  }
  dat=nbfit_tree_setup_stan_data(tree,agedf,rate_switch_nodes,xcross,b_pool_rates =b_pool_rates)
  dat$early_growth_model_on=early_growth_model_on
  dat$concentration=rep(concentration,dat$NINT)
  dat$q=dat$q[dat$xidx>0]
  stanr=rstan::sampling(stanmodels[[model]], data=dat,chains = cores,iter=niter,cores=cores,control=stan_control)
  if(length(split_nodes)>0){
    post.dist=rstan::extract(stanr)
    idx.rate.switch=match(split_nodes,tree$edge[,2])
    tmp=tree
    nodetdist=matrix(apply(post.dist$ta,1,function(x){ tmp$edge.length=x;nodeHeights(tmp)[idx.rate.switch,2]}),nrow=length(idx.rate.switch))
    nodetdistl=matrix(apply(post.dist$ta,1,function(x){ tmp$edge.length=x;nodeHeights(tmp)[idx.rate.switch,1]}),nrow=length(idx.rate.switch))
    ttmp=apply(nodetdist,1,function(x) quantile(x,prob=c(0.025,0.5,0.975),na.rm=TRUE))
    nodet=list(lb95=ttmp[1,],ub95=ttmp[3,],median=ttmp[2,],mean=apply(nodetdist,1,mean))
    ttmp=apply(nodetdistl,1,function(x) quantile(x,prob=c(0.025,0.5,0.975,na.rm=TRUE)))
    nodetl=list(lb95=ttmp[1,],ub95=ttmp[3,],median=ttmp[2,],mean=apply(nodetdistl,1,mean))
  }else{
    nodet=NULL
    nodetl=NULL
  }
  sumc=rstan::summary(stanr)$summary
  idx=grepl("^lambda",rownames(sumc))
  lambda=list(mean=sumc[idx,"mean"],sd=sumc[idx,"sd"],lb=sumc[idx,"2.5%"],ub=sumc[idx,"97.5%"],median=sumc[idx,"50%"])
  idx=grepl("^p",rownames(sumc))
  p=list(mean=sumc[idx,"mean"],sd=sumc[idx,"sd"],lb=sumc[idx,"2.5%"],ub=sumc[idx,"97.5%"],median=sumc[idx,"50%"])
  idx=grepl("^k",rownames(sumc))
  k=list(mean=sumc[idx,"mean"],sd=sumc[idx,"sd"],lb=sumc[idx,"2.5%"],ub=sumc[idx,"97.5%"],median=sumc[idx,"50%"])
  idx=grepl("^ta",rownames(sumc))
  outtree=tree
  outtree$edge.length=sumc[idx,"mean"]
  list(fullres=stanr,lambda=lambda,p=p,k=k,ultratree=outtree,intree=tree,dat=dat,upper_node_lims=nodet,lower_node_lims=nodetl,split_nodes=split_nodes)
}

#' Sets up input data for the stan models.
nbfit_tree_setup_stan_data=function(tree,##<< tree with mutation counts per branch. Best to collapse polytomies..
         tip_heights,##<< dataframe with tip.label, age, specifying the terminal timepoint for each tip.
         rate_switch_nodes,##<< child node ids of edges with a switch in lambda.
         xcross=rep(0.95,length(rate_switch_nodes)),##<< How far down the branch the switch occurs.
         b_pool_rates=FALSE
){
  if(is.null(tip_heights$tip.label)){
    stop("tip_heights parameter should contain an age and tip.label field")
  }
  L=length(tree$edge.length)
  N=length(tree$tip.label)
  ##The following gives the index of the tips in the edge matrix
  tip_heights=tip_heights[match(tree$tip.label,tip_heights$tip.label),]
  ##parallel to edge. This gives the minimum sampling age of tips descendent from specified node.
  tip_min_age=sapply(tree$edge[,2],function(x) min(tip_heights$age[intersect(getDescendants(tree,x),1:N)]))
  ###
  parentidx=match(tree$edge[,1],tree$edge[,2])-1
  parentidx=ifelse(is.na(parentidx),-1,parentidx)
  #Index of tips in edge matrix
  tipidx=match(1:N,tree$edge[,2])#(tree$edge[,2]<=N)
  #Index of interior edges in edges matrix # This is the key for entities that live in NINT space.
  xidx=setdiff(1:L,tipidx)
  ##map edges to the interior edge matrix.
  pxidx=match(1:L,xidx)
  #map node ids to edge indices...  Root will have NA index
  pidx=match(1:(L+1),tree$edge[,2])
  LX=length(xidx)##L-N
  edges=tree$edge
  #mutation count per edge
  m=tree$edge.length
  s=tree$sensitivity
  ##Get ball park lambda
  tree2=tree
  tree2$edge.length=tree2$edge.length/tree$sensitivity
  ##Adjusted node heights
  nh=nodeHeights(tree2)

  tree3=tree
  tree3$edge.length=rep(1,length(tree2$edge.length))
  ##Node heights as number of ancestors
  mh=nodeHeights(tree3)

  idx.tip=match(1:length(tree$tip.label),tree$edge[,2])
  lambda_est=median(nh[idx.tip,2]/tip_min_age[idx.tip])
  cat(sprintf("Median lambda estimate=%3.2f\n",lambda_est))
  if(is.na(lambda_est)){
    stop("Missing lambda estimate!")
  }
  ##Parents of every edge including the edge itself [order from root to tip]
  parentlist=lapply(edges[,2],function(node) {
    ##We need to make sure that parent edges are ordered from root to tip
    ##TODO export get_parents in rsimpop
    parents=rev(rsimpop:::get_parents(node,edges))##
    if(length(parents)>0){
      #parents
      match(parents,edges[,2])
    }else{
      parents
    }
  }
  )
  ##Index of applicable lambda rate - 0 for wild type.
  rates=rep(0,length(m))
  if(length(rate_switch_nodes)>0){
    for(i in 1:length(rate_switch_nodes)){
      snode=rate_switch_nodes[i]
      nodes=c(snode,get_all_node_children(node = snode,tree = tree))
      ##Check that nodes in correct order.
      if(length(unique(rates[match(nodes,tree$edge[,2])]))>1){
        stop("switch nodes provided in wrong order!")
      }
      idxanc=match(nodes,tree$edge[,2])
      if(length(idxanc)>0){
        if(b_pool_rates){
          ##do we just want to estimate one extra rate?
          rates[idxanc]=1
        }else{
          rates[idxanc]=i
        }
      }
    }
  }
  #Rates of parent node.
  ratesp=rates[match(tree$edge[,1],tree$edge[,2])]
  #If node has no parent assume the wild type rate
  ratesp[which(is.na(ratesp))]=0
  #Index of edges where rate switches occur.
  nlambda=max(rates)+1
  #alpha=rep(0,L)
  #if(length(rate_switch_nodes)>0){
  #  alpha[match(rate_switch_nodes,tree$edge[,2])]=xcross
  #}
  ##Get approximately ultrametric tree
  aut=make_almost_ultrametric(tree2,tip_heights)
  ## prior estimate for t..
  tprior=aut$edge.length/tip_min_age
  parentidx=parentidx+1
  q=rep(NA,L)
  ##shape2v=rep(-1,L)
  for(i in 1:L){
    k=i;##parentidx[i];
    ptot=0
    while(parentidx[k]>0){
      k=parentidx[k];
      ptot=ptot+tprior[k]
    }
    q[i]=tprior[i]/(1-ptot)
  }
  if(length(rate_switch_nodes)<=1){
    idxcrossover=as.array(match(rate_switch_nodes,tree$edge[,2]))
  }else{
    idxcrossover=match(rate_switch_nodes,tree$edge[,2])
  }
  list(N=L,#Number of branches
       NINT=LX,#Number of internal branches
       NLAMBDA=nlambda,#Number of rates being inferred (includes wild type)
       parentidx=parentidx,#Index of parent branch (-1 for root)
       xidx=ifelse(is.na(pxidx),0,pxidx),#Index in internal edge matrix (0 if not internal)
       tip_min_age=tip_min_age,# For each branch this gives the minimum sampling age if descendent branches
       rates=rates+1,# Lambda index of rate at end of branch i.e. in range 1:NLAMBDA
       ratesp=ratesp+1,# Lambda index of rate of parent branch
       m=m,# Observed mutation count on branch
       s=s,# Estimated sensitivity of mutation detection on branch.
       lambda_est=lambda_est,# Initial estimate of lambda (used to establish a weak prior)
       nh=mh[xidx[xidx>0],2],# Heights of internal nodes [TODO get rid of redundant xidx>0 criteria ]
       q=q,# Prior for stick breaking fraction.
       idxcrossover=idxcrossover # index of branches where the lambda change point occurs (excludes wild type)
  )
}

##Approximate approach making a tree almost timepoint specific
##ultrametric and guaranteeing that root to tip is less than or equal to sample age
make_almost_ultrametric=function(tree,agedf){
  ## min age  of sampling
  tma=sapply(tree$edge[,2],
             function(x) min(agedf$age[intersect(getDescendants(tree,x),1:length(tree$tip.label))]))
  tree$edge.length=tree$edge.length+1 ## Add one
  nh=nodeHeights(tree)  ## Cumulative root to branch mutation count
  ## Maximum mutation count of tips that are descendent from each branch
  nha=sapply(tree$edge[,2],function(x) max(nh[match(getDescendants(tree,x),tree$edge[,2]),2]))
  ut=tree
  ut$edge.length=(tree$edge.length/nha)*tma
  ut
}


