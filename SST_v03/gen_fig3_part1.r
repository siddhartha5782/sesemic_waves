#
# Run TBATS on Y0 (two harmonic components)
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
    if (!is.null(x$seasonal.periods)){
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

	yhat <- yhat + rowSums(out[,3:ncol(out)])
    }	
    else{
	len = length(yhat);
	seasonality = 0*sign(1:len);
	out <- cbind(out, seasonality); 
	out <- cbind(out, seasonality); 

        colnames(out)[3:2] <- paste("season",1:2,sep="")
    }

    residuals <- y - yhat
    out <- cbind(out, yhat=yhat, residuals=residuals, errors=x$errors)
	
# Add time series characteristics
    out <- ts(out)
    tsp(out) <- tsp(y)
	
# Do the plot
    out
}



###############################################
    # number of realizations
realizationNo = 51;
simulation = 0;	

rm(x)

n = 1000;
t = (1:n)/100;


    # for the first simulation
s1 = 2.5*cos(2*pi*t);
s2 = 3*cos(2*pi*pi*t);

s = s1 + s2;

trend = 8*(1/(1+(t/5)^2)+exp(-t/10));

    # noise level. K = 2 or 5 in simulation 1
dt = 1+0.1*cos(pi*t);

    # the L2 error of the estimator
err_s1 = 1:realizationNo;
err_s2 = 1:realizationNo;
err_arma = 1:realizationNo;
err_trend = 1:realizationNo;
runtime = 1:realizationNo;
recon_s1 = 1:n;
recon_s2 = 1:n;
recon_arma = 1:n;
recon_trend = 1:n;
estPeriods = matrix(nrow=1, ncol=2);
allx = 1:n;



for (ii in 1:realizationNo) { 

    set.seed(12101+ii)

    if(simulation==1){
	    # locally stationary noise
        Yt1 = arima.sim(list(order=c(1,0,1), ar=c(-0.5), ma=c(0.4)), n=n/2, innov=rt(n/2,4));
        Yt2 = arima.sim(list(order=c(1,0,1), ar=c(0.2), ma=c(0.51)), n=n/2, innov=rt(n/2,4));
        Yt = c(4*Yt1,Yt2)
        noise = dt*Yt; 
    }

    if(simulation==0){
	    # previous simulation. 
        Yt = arima.sim(list(order=c(1,0,1), ar=c(-0.5), ma=c(0.4)), n=n, innov=rt(n,4));
        noise = 2*dt*Yt;
    } 

    x = s + trend + noise;
    allx = cbind(allx, x);


        # run TBATS
    ptm = proc.time();

    xshift = -min(x)+1;
    AICvalue = matrix(nrow=s1valNo, ncol=s2valNo);

    tmp_p1 = 100;
    tmp_p2 = 100/pi;

    x.tbats <- tbats(x+xshift, use.box.cox=FALSE, use.trend=TRUE, use.damped.trend=NULL, seasonal.periods=c(100, 100/pi), use.arma.errors=FALSE);

    tmp = proc.time() - ptm;
    runtime[ii] = tmp[1];

 
        # extract seasonality and trend
    x.tbats.out <- out.tbats(x.tbats);


	# the error of the i-th realization
    err_s1[ii] = sqrt(sum((x.tbats.out[,4] - s1)^2))/sqrt(sum(s1^2));
    err_s2[ii] = sqrt(sum((x.tbats.out[,3] - s2)^2))/sqrt(sum(s2^2));
    err_arma[ii] = sqrt(sum((x.tbats.out[,6] - noise)^2))/sqrt(sum(noise^2));
    err_trend[ii] = sqrt(sum((x.tbats.out[,2] - xshift - trend)^2))/sqrt(sum(trend^2));

	
	# record all the estimated seasonalities
    recon_s1 = cbind(recon_s1, x.tbats.out[,4]);
    recon_s2 = cbind(recon_s2, x.tbats.out[,3]);
    recon_arma = cbind(recon_arma, x.tbats.out[,6]);
    recon_trend = cbind(recon_trend, x.tbats.out[,2] - xshift);
   

	# record all the estimated seasonal periods
    tmpPeriod = matrix(nrow=1, ncol=2);
    if(!is.null(x.tbats$seasonal.periods)){
        tmpPeriod[1] = 100/(s1start + s1jump*optjj);
        tmpPeriod[2] = 100/(s2start + s2jump*optkk);
        estPeriods = rbind(estPeriods, tmpPeriod);
    }else{
        tmpPeriod[1] = -1;
        tmpPeriod[2] = -1;
        estPeriods = rbind(estPeriods, tmpPeriod);
    }
}


    # the realization with the   mediai error
OptIdx = which(err_s1 == median(err_s1), arr.ind = TRUE);
final_s1 = recon_s1[,OptIdx+1];
final_s2 = recon_s2[,OptIdx+1];
final_arma = recon_arma[,OptIdx+1]; 
final_trend = recon_trend[,OptIdx+1];

    # plot the results
dev.new()
par(mfrow = c(4,1), mar=c(4,4,1,1)+0.1)

plot(t, x, type="l", xaxt="n", xlab="", ylab=expression(X(t)));

plot(t, final_trend, type="l", lty=1, col="black", xaxt="n", xlab="", ylab="");
par(new=TRUE);
plot(t, trend, type="l", pch=22, lty=2, col=rgb(r=0.5,g=0.5,b=0.5));

plot(t, final_s1, type="l", lty=1, col="black", xaxt="n", xlab="", ylab="");
lines(t, s1, type="l", pch=22, lty=2, col=rgb(r=0.5,g=0.5,b=0.5));

plot(t, final_s2, type="l", lty=1, col="black", xaxt="n", xlab="", ylab="");
lines(t, s2, type="l", pch=22, lty=2, col=rgb(r=0.5,g=0.5,b=0.5));

report = cbind(s1Mean = mean(err_s1), s1SD = sd(err_s1), s2Mean = mean(err_s2), s2SD = sd(err_s2), trendMean = mean(err_trend), trendSD = sd(err_trend), armaMean = mean(err_arma), armaSD = sd(err_arma), timeMean = mean(runtime), timeSD = sd(runtime));

x1hat = final_s1;
x2hat = final_s2;
trendhat = final_trend;
noisehat = final_arma;


xt = allx[,OptIdx+1];
x1 = s1;
x2 = s2;

if(simulation==1)
    exportdata1 = data.frame(xt,t,x1,x2,trend,noise,x1hat,x2hat,trendhat,noisehat);

if(simulation==0)
    QQQ1 = data.frame(xt,t,x1,x2,trend,noise,x1hat,x2hat,trendhat,noisehat);

