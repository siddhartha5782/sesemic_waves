require(astsa)
require(forecast)



    # sr = sampling rate
sr = 100;

    # how many points
n = 1000

    # time
time = (1:n)/sr;

    # the A_{\epsilon,d} class function
am1 = (1+0.1*cos(time))*rev((atan((1:n)/43.5-10)))/2 + 2;
am2 = 3.5*sign(1:n); am2[751:1000] = 2*sign(1:250);
phase1 = time + 0.1*sin(time);
if1 = 1 + 0.1*cos(time);
phase2 = 3.4*time - 0.02*time^(2.3);
if2 = 3.4 - 0.069*time^(1.3);

s21 = am1*cos(2*pi*phase1)
s22 = am2*cos(2*pi*phase2);
s2 = s21+s22;

trend = 8*(1/(1+(time/5)^2) + exp(-time/10));
trend2 = 2*time + 10*exp(-((time-4)^2)/6);


dt = 1 + 0.1*cos(pi*time);
K = 5;

clean = s2 + trend;
clean2 = s2 + trend2;
noise = arima.sim(list(order=c(1,0,1), ar=c(-0.5), ma=c(0.4)), n=n, sd=K)
x = clean + dt*noise;


par(mfrow = c(5,2), mar=c(0.1,5,0.1,0.5)+0.1, oma=c(2,1,1,1), cex=0.8);
plot(time, am1, type="l", lty=1, col="black", xaxt="n", xlab="", ylab=expression(A[1]),ylim=c(0,4));
plot(time, am2, type="l", lty=1, col="black", xaxt="n", xlab="", ylab=expression(A[2]),ylim=c(0,4));
plot(time, if1, type="l", lty=1, col="black", xaxt="n", xlab="", ylab=expression(phi[1]^"'"),ylim=c(0,4));
plot(time, if2, type="l", lty=1, col="black", xaxt="n", xlab="", ylab=expression(phi[2]^"'"),ylim=c(0,4));
plot(time, s21, type="l", lty=1, col="black", xaxt="n", xlab="", ylab=expression(s[21]));
plot(time, s22, type="l", lty=1, col="black", xaxt="n", xlab="", ylab=expression(s[22]));
plot(time, trend, type="l", lty=1, col="black", xaxt="n", xlab="", ylab=expression(T[1]),ylim=c(0,20));
plot(time, clean, type="l", lty=1, col="black", xaxt="n", xlab="", ylab=expression(s[2]+T[1]),ylim=c(-5,25));
plot(time, trend2, type="l", lty=1, col="black", xlab="Time", ylab=expression(T[2]),ylim=c(0,20));
plot(time, clean2, type="l", lty=1, col="black", xlab="Time", ylab=expression(s[2]+T[2]),ylim=c(-5,25));

