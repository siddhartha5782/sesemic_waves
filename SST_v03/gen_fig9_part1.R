#
# run TBATS on varicella and HZ. **The results are not saved** 
#
require(astsa)
require(forecast)


out.tbats <- function (x, main="Decomposition by TBATS model", ...) 
{
# modified from plot.tbats in tbats.R	
# Get original data, transform if necessary
    if (!is.null(x$lambda)) 
    	y <- InvBoxCox(x$y, x$lambda)
    else 
    	y <- x$y

# Compute matrices
    tau <- ifelse(!is.null(x$k.vector), 2*sum(x$k.vector), 0)
    w <- .Call("makeTBATSWMatrix", smallPhi_s = x$damping.parameter, kVector_s=as.integer(x$k.vector), arCoefs_s = x$ar.coefficients, maCoefs_s = x$ma.coefficients, tau_s=as.integer(tau), PACKAGE = "forecast")
	
    out <- cbind(observed=c(y), level=x$x[1,])
    nonseas <- 1+!is.null(x$beta) # No. non-seasonal columns in out

# out[,1] observation;  out[,2] level;  out[,3:4] seasonality
    yhat <- out[,2]
	
# number of seasonality
    nseas <- length(x$seasonal.periods) # No. seasonal periods

    seas.states <- cbind(x$seed.states,x$x)[-(1:(1+!is.null(x$beta))),]
    seas.states <- seas.states[,-ncol(seas.states)]
    w <- w$w.transpose[,-(1:(1+!is.null(x$beta))),drop=FALSE]
    w <- w[,1:tau,drop=FALSE]
    j <- cumsum(c(1,2*x$k.vector))

    for(i in 1:nseas)
        out <- cbind(out, season=c(w[,j[i]:(j[i+1]-1),drop=FALSE] %*% seas.states[j[i]:(j[i+1]-1),]))

    if(nseas > 1)
        colnames(out)[nonseas + 1:nseas] <- paste("season",1:nseas,sep="")

    if(nseas > 1)
        yhat <- yhat + rowSums(out[,3:ncol(out)])
    else
        yhat <- yhat + out[,3]
	

    residuals <- y - yhat
    out <- cbind(out, yhat=yhat, residuals=residuals, errors=x$errors)
	
# Add time series characteristics
    out <- ts(out)
    tsp(out) <- tsp(y)
	
# Do the plot
    #plot(out, main=main, nc=1, ...)
    out
}


dirnm <- "/Users/hauwu/Dropbox/tex20120723/_WORKING_seasonality/fig/"
filenm <- paste(dirnm, "Varicella.20120127.1mCohort.txt", sep="")
varicella <- matrix(scan(filenm), byrow=T,ncol=4)
varicella <- varicella[,4] / varicella[,2]


###############################################
    # analyze varicella
n = 520;
t = (1:n)/52+2000;

x = varicella*1000
xshift = -min(x)+1;


x.tbats <- tbats(x+xshift, use.box.cox=FALSE, use.trend=TRUE, use.damped.trend=NULL, seasonal.periods=c(52), use.arma.errors=FALSE);

    # extract seasonality and trend
x.tbats.out <- out.tbats(x.tbats);

s1V = x.tbats.out[,3];
trendV = x.tbats.out[,2] - xshift;
armaV = varicella*1000 - s1V - trendV;


###############################################
    # analyze HZ

filenm <- paste(dirnm, "500.HZ.20120127.1mCohort.txt", sep="")
HZ <- matrix(scan(filenm), byrow=T,ncol=4)
HZ <- HZ[,4] / HZ[,2]

n = 520;
t = (1:n)/52+2000;

y = HZ*10000
yshift = -min(y)+1;

y.tbats <- tbats(y+yshift, use.box.cox=FALSE, use.trend=TRUE, use.damped.trend=NULL, seasonal.periods=c(52), use.arma.errors=FALSE);

    # extract seasonality and trend
y.tbats.out <- out.tbats(y.tbats);

s1H = y.tbats.out[,3];
trendH = y.tbats.out[,2] - yshift;
armaH = HZ*10000 - s1H - trendH;

