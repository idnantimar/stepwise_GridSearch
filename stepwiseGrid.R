library(NMOF)
library(corrplot)
library(parallel)
library(xfun)


gridSearch_Parallel<-function(MINIMIZE,currentGrid,clust=makeCluster(detectCores()-1),export_data=NULL,needed_packages=list("base")){

  clusterExport(clust,export_data)
  clusterExport(clust,"needed_packages",envir = environment())
  clusterEvalQ(clust,xfun::pkg_attach(needed_packages))
  
  expandGrid = expand.grid(currentGrid)
  VALUES = parApply(clust, X = expandGrid, MARGIN = 1,FUN = MINIMIZE) ## the main computation chunk
  optVal = min(VALUES)
  opt_ix = which(VALUES==optVal)
  if(length(opt_ix)>1) opt_ix = sample(opt_ix,1) # if multiple minima, randomly choose one
  
  Out = list()
  Out[['minfun']] = optVal
  Out[['minlevels']] = as.numeric(expandGrid[opt_ix,])
  Out[['values']] = VALUES
  return(Out)
}




stepwiseGrid<-function(FUN,
                initialGrid = list(seq(0,1,0.2),seq(0,1,0.2)),
                minimize = TRUE ,
                updateMesh=2,update=1,
                plotting=FALSE,
                para = list(detectCores()-1, ls(envir = .GlobalEnv),loadedNamespaces()) ){
  
  
  ## FUN : the objective function , to be minimized over the grid of hyperparameters
  ## initialGrid : possible choices of hyperparameters to start with
              #  a list where each element is a sequence corresponding to a hyperpameter
  ## minimize : whether it is a minimization(default) or maximization(input FALSE) problem
  ## updateMesh : a single number k or a vector (k1,k2,...)
              # new k1 points will be added between the current optimum and the midway to adjacent choice along 1st hyperparameter, k2 along 2nd one, ...
              # so, updated grid is of size (2k1+3)*(2k2+3)*...
  ## update : number of times the grid to be updated
          # put 0 if initial grid is enough
  ## plotting : possible function values will be plotted (upto 2D case only)
  ## para : whether to use parallel computation or not
        # input is a list containing - 'number of clusters' & 'list of variables and functions to be exported' , 'list of libraries in use as list("package_name1","package_name2",...)'
        # If not sure about input, don't change the default one. But using appropriate input will be efficient
        # to disable parallel computation input NULL


MINIMIZE = ifelse(minimize,FUN, function(X) -1*FUN(X))
if(is.null(para)){ 
  METHOD = function(MINIMIZE,currentGrid) NMOF::gridSearch(MINIMIZE,currentGrid,method = 'loop',printDetail = F)
}else{
  clust=makeCluster(para[[1]])
  export = para[[2]]
  packages = para[[3]]
  METHOD = function(MINIMIZE,currentGrid) gridSearch_Parallel(MINIMIZE,currentGrid,clust,export,packages)
} 


soln = METHOD(MINIMIZE,initialGrid)
optimum = soln$minfun
optPar = soln$minlevels
dimPar = length(optPar)


if(plotting*(dimPar==2)){
  plot_it<-function(currentGrid,currentSoln,step){
    cat(c("\n",paste("Step",step),rep("..",step+1),"\n"),sep = "")
    
    try({r = length(currentGrid[[1]])
    c = length(currentGrid[[2]])
    val = matrix(currentSoln$values,nrow=r)*ifelse(minimize,1,-1)
    rownames(val) = currentGrid[[1]]
    colnames(val) = currentGrid[[2]]
    corrplot(val,is.corr = F,method = 'color')
    mtext(paste('step',step,ifelse(minimize," | min"," | max")),side=ifelse(r>=c,4,1))})
  }
  plot_it(initialGrid,soln,0)
  
}else if(plotting*(dimPar==1)){
  plot_it<-function(currentGrid,currentSoln,step){
    cat(c("\n",paste("Step",step),rep("..",step+1),"\n"),sep = "")
    
     try({val = ifelse(minimize,1,-1)*currentSoln$values
    plot(currentGrid[[1]], val,type = 'b',col='blue',xlab='hyperparameter',ylab = 'value',main=paste('step',step,ifelse(minimize," | min"," | max")))})
  }
  plot_it(initialGrid,soln,0)
  
}else{
  plot_it<-function(currentGrid,currentSoln,step){
    cat(c("\n",paste("Step",step),rep("..",step+1),"\n"),sep = "")
   }
  plot_it(initialGrid,soln,0)
  
} 
DATA = list()
stepData<-function(currentSoln){
  val = list()
  val[["value"]] = ifelse(minimize,1,-1)*currentSoln$minfun
  val[["select"]] = currentSoln$minlevels

  return(val)
}
DATA[[paste("step ",0)]] = stepData(soln)


step = 1
currentGrid = initialGrid
K = 2+ rep(1,dimPar)*updateMesh 
  # in between the midway of optimum choice and adjacent choice , K-2 new possible choices will be added at every direction 
    # current choice is also there in the new grid 
while(step<=update){
  currentPar = optPar
  ixOpt = sapply(1:dimPar,function(t) which(currentPar[t]== currentGrid[[t]]))
  currentGrid = lapply(1:dimPar, 
                       function(t){ 
                         A = currentGrid[[t]]
                         ix = ixOpt[t]
                         k = K[t]
                         return(union(rev(seq(A[ix],(A[ix]+A[max(ix-1,1)])/2,length.out=k)),
                               seq(A[ix],(A[min(ix+1,length(A))]+A[ix])/2,length.out=k)))
                         })  
  
  soln = METHOD(MINIMIZE,currentGrid)
  optimum = soln$minfun
  optPar = soln$minlevels
    
  plot_it(currentGrid,soln,step)
  DATA[[paste("step ",step)]] = stepData(soln)
    
  step = step + 1
  
}


if(!is.null(para)) stopCluster(clust)

OUT = list()
OUT[['optimumHypereparameter']] = optPar
OUT[[ifelse(minimize,'minValue','maxValue')]] = optimum*ifelse(minimize,1,-1)
OUT[["stepwise_Data"]] = DATA
return(OUT)
  ## returns - 'optimum choice of hyperparameters' & 'corresponding optimized value of the function' 

}
