
############################# BIENVENUE #############################

############# I SOLEMNLY SWEAR THAT I AM UP TO NO GOOD ##############

#################################
## Loading packages

library(MASS)
library(lavaan)

################################
## The NCArigour function
## Estimates expected necessity effect due to correlation between X and Y
## The following values are plugged into the function:
## n = sample size
## rxy = empirical correlation between X and Y
## dnecc = empirical necessity effect
## ci.req = requested confidence interval (default = 0.95)
## cr_fdh = T if estimating CR-FDH (default) and cr_fdh = F if estimating CE-FDH
## R = number of simulated datasets (default = 1000)
## silent = F if printing results and silent = T if not

NCArigour <- function(n,rxy,dnecc,ci.req=0.95,cr_fdh=T,R=1000,silent=F){
  
  rm <- matrix(c(1,rxy,rxy,1),nrow=2) ## correlation matrix
  
  simres <- unlist(sapply(1:R, function(i){ ## for each simulation
    set.seed(i) ## for reproducibility
    df <- data.frame(mvrnorm(n=n,mu=c(0,0),Sigma=rm,empirical=T)) ## the data
    
    names(df) <- c("x","y") ## naming variables in data
    
    x.max <- max(df$x)
    y.max <- max(df$y)
    x.min <- min(df$x)
    y.min <- min(df$y)
    
    scope <- (x.max-x.min)*(y.max-y.min)
    
    sx <- df$x[order(df$x)] ## sorted X
    sy <- df$y[order(df$x)] ## Y in the order of sorted X
    
    cmy <- cummax(sy) ## cumulative max on the Y-variable
    cp <- ifelse(sy>=cmy,1,0) ## if the Y-value is a ceiling point
    
    cx <- sx[cp==1] ## X-values for ceiling points
    cy <- sy[cp==1] ## Y-values for ceiling points
    
    minx <- which(df$x==min(df$x))
    maxy <- which(df$y==max(df$y))
    same <- sum(minx %in% maxy)
    
    if(same>0) nec <- 0
    
    if(same==0){
      
      if(cr_fdh==F){
        
        numc <- length(cx)
        x.base <- cx[2:numc] - cx[1:(numc-1)]
        y.height <- y.max - cy[1:(numc-1)]
        nec <- sum(x.base*y.height) / scope
      }
      
      if(cr_fdh==T){
        
        fitc <- lm(cy ~ cx)
        ylim <- predict(fitc,newdata=data.frame(cx=x.min))
        xlim <- (y.max-fitc$coefficients[1])/fitc$coefficients[2]
        nec <- (0.5 * (y.max-ylim) * (xlim-x.min)) / scope
      }
    }
    
    return(nec) ## returns the necessity effect
  }))
  
  d <- mean(simres, na.rm=T) ## mean/expected necessity effect across the R effects
  p1 <- (1-ci.req)/2 ## lower percentile
  p2 <- 1-p1 ## upper percentile
  ci <- quantile(simres, probs=c(p1,p2)) ## required percentile confidence interval
  pv <- round(sum(simres>dnecc, na.rm=T)/R,3) ## frequency of effects stronger than the empirical effect
  tx1 <- ifelse(pv<0.05,"null rejected","null retained") ## text if p < 0.05 or not
  tx2 <- ifelse(ci[2]<dnecc,"null rejected","null retained") ## text if empirical effect above range of spur. or not
  bm <- round(max(simres,na.rm=T),3) ## max spurious effect
  
  if(silent==F){ ## if printing the results
    
    cat("\n\n*****",
        "\nd(obs)=",round(dnecc,3),
        "\nmax(sim)=",bm,
        "\np(sim>d(obs))=",pv,
        "\n",tx1,
        "\n**",
        "\nd(sim)=",round(d,3),
        "\nci(percentile)=[",round(ci[1],3),",",round(ci[2],3),"]",
        "\n",tx2,
        "\n*****\n")
  }
  
  res <- c(d,ci,pv) ## what the function returns
  return(res)
}

#################################
## Data

n <- 326 ## sample size in Yang et al. (2026)

orvar <- c("OS","PSCul","AERP","PI","SE","PSC") ## order of variables

rm <- matrix(c(
  
  1.000, 0.463, 0.278, 0.380, 0.430, 0.482,
  0.463, 1.000, 0.316, 0.434, 0.482, 0.563,
  0.278, 0.316, 1.000, 0.526, 0.401, 0.484,
  0.380, 0.434, 0.526, 1.000, 0.405, 0.576,
  0.430, 0.482, 0.401, 0.405, 1.000, 0.637,
  0.482, 0.563, 0.484, 0.576, 0.637, 1.000),nrow=6)

rownames(rm) <- colnames(rm) <- orvar

#################################
## Figure 2

## Order of variables: orvar (see above)

rval <- rm[6,1:5] ## correlations with PSC

yl <- rev(orvar[1:5]) ## y-labels

numef <- length(yl) ## number of effects
upplim <- numef+0.5 ## adding some space

panlab <- c("A: CE-FDH", "B: CR-FDH") ## panel labels

sf <- 0.8 ## sizing factor

if(dev.cur()==2) dev.off() ## removing earlier plots
par(mar=c(1,1,1.2,0), oma=c(2.5,2.5,0.5,1), mfrow=c(1,2))

for (j in 1:2){
  
  if(j==1){ ## first CE-FDH
    dval <- c(0.221,0.232,0.105,0.246,0.236) ## effects 
    crval <- rep(0,5)
  }
  
  if(j==2){ ## then CR-FDH
    dval <- c(0.197,0.216,0.108,0.242,0.218) ## effects
    crval <- rep(1,5)
  }
  
  plot(c(0,0.5), c(0.5,upplim), type="n", xlab="", ylab="", yaxt="n",cex.axis=sf, las=1) ## empty plot
  axis(2,1:numef, labels=F,las=2, cex.axis=sf) ## Y-axis
  if(j==1) axis(2,1:numef, labels=yl,las=2, cex.axis=sf) ## Y-axis
  
  if(j==1) mtext("Necessity Effect", side = 1, line=1.2, outer=T, cex=sf) ## adding X-label
  mtext(panlab[j],3, line=0.2, cex=sf) ## panel label
  
  resemp <- t(sapply(1:numef,function(i){ ## for each necessity effect
    
    res <- NCArigour(n=n,rxy=rval[i],dnecc=dval[i],cr_fdh=crval[i],silent=F) ## applying the NCArigour function
    
    yval <- (numef+1)-i ## Y-value (reversed to start on top)
    
    abline(h=yval, col="gray") ## a thin gray line
    
    arrows(res[2],yval,res[3],yval,angle=90,length=0.05, code=3,lwd=2) ## range of spuriousness
    points(dval[i], yval, pch=21, col="black",bg="black") ## the empirical necessity effect
    
    return(res) ## returning the results by the NCArigour function
    
  }))
}

#################################
## Alternative model

spurmod <- "

## Loadings

GSE =~ OS+PSCul+AERP+PI+SE+PSC

## Correlated residuals

AERP ~~ PI
PI ~~ SE
OS ~~ PSCul

## (Error) variances

GSE ~~ 1*GSE
OS ~~ OS
PSCul ~~ PSCul
AERP ~~ AERP
PI ~~ PI
SE ~~ SE
PSC ~~ PSC

"

fit.spur <- lavaan(spurmod, sample.cov=rm, ## fitting model to data 
              sample.nobs=n, meanstructure=F)

summary(fit.spur, fit.measures=T, ci=T, standardized=T, rsq=T) ## the results


########################## MISCHIEF MANAGED #########################

############################# AU REVOIR #############################

