require(astsa)
require(forecast)
require(R.matlab)


    # number of realizations
realizationNo = 201;

n = 1000;
t = (1:n)/100;

dt = 1+0.1*cos(pi*t);
allx0 = 1:n;
allx1 = 1:n;
allx1g = 1:n;
allx1c = 1:n;
allx1tg = 1:n;

set.seed(12301)

for (ii in 1:realizationNo) { 
    Yt1 = arima.sim(list(order=c(1,0,1), ar=c(-0.5), ma=c(0.4)), n=n/2, innov=rt(n/2,4));
    Yt2 = arima.sim(list(order=c(1,0,1), ar=c(0.2), ma=c(0.51)), n=n/2, innov=rt(n/2,4));
    Yt = c(5*Yt1,Yt2)
    noise = dt*Yt;
    allx1 = cbind(allx1, noise);

    Yt1g = arima.sim(list(order=c(1,0,1), ar=c(-0.5), ma=c(0.4)), n=n/2, innov=rnorm(n/2));
    Yt2g = arima.sim(list(order=c(1,0,1), ar=c(0.2), ma=c(0.51)), n=n/2, innov=rnorm(n/2));
    Ytg = c(5*Yt1g,Yt2g)
    noiseg = dt*Ytg;
    allx1g = cbind(allx1g, noiseg);

    Yttg = c(5*Yt1,Yt2g)
    noisetg = dt*Yttg;
    allx1tg = cbind(allx1tg, noisetg);



    Yt1c = arima.sim(list(order=c(1,0,1), ar=c(-0.5), ma=c(0.4)), n=n/2, innov=rt(n/2, 1));
    Yt2c = arima.sim(list(order=c(1,0,1), ar=c(0.2), ma=c(0.51)), n=n/2, innov=rt(n/2, 1));
    Ytc = c(5*Yt1c,Yt2c)
    noisec = dt*Ytc;
    allx1c = cbind(allx1c, noisec);


    Yt0 = arima.sim(list(order=c(1,0,1), ar=c(-0.5), ma=c(0.4)), n=n, innov=rt(n,4));
    noise = 2*dt*Yt0; 
    allx0 = cbind(allx0, noise);
}

writeMat("t4_sim0_noise.mat", Rdata = allx0);
writeMat("t4_sim1_noise.mat", Rdata = allx1);
writeMat("gaussian_sim1_noise.mat", Rdata = allx1g);
writeMat("cauchy_sim1_noise.mat", Rdata = allx1c);
writeMat("tg_sim1_noise.mat", Rdata = allx1tg);

