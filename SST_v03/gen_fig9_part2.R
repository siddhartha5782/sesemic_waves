#
# display the TBATS results of varicella and HZ. 
#


###########################
    # plot the results
#par(mfrow = c(4,2), mar=c(0.1,1.5,0.8,0.5)+0.1, oma=c(1,1,1,2), cex=0.9);
par(mfrow = c(4,2));
par(mar=c(0.1,2,1,0.8)+0.1, oma=c(2,1,1,1), cex=1.2)


xrange = range(varicella*1000);
yrange = range(HZ*10000);

plot(t, varicella*1000, type="l", xaxt="n", xlab="", ylab="",ylim=xrange);
par(new=TRUE); plot(t, trendV+s1V, type="l", col="green", xaxt="n", yaxt="n", xlab="", ylab="" ,ylim=xrange)

plot(t, HZ*10000, type="l", xaxt="n", xlab="", ylab="",ylim=yrange);
par(new=TRUE); plot(t, trendH+s1H, type="l", col="green", xaxt="n", yaxt="n", xlab="", ylab="" ,ylim=yrange)

plot(t, varicella*1000, type="l", lty=1, col="black", xaxt="n", xlab="", ylab="",ylim=xrange);
par(new=TRUE); plot(t, trendV, type="l", col="red", xaxt="n", yaxt="n", xlab="", ylab="" ,ylim=xrange)

plot(t, HZ*10000, type="l", lty=1, col="black", xaxt="n", xlab="", ylab="",ylim=yrange);
par(new=TRUE); plot(t, trendH, type="l", col="red", xaxt="n", yaxt="n", xlab="", ylab="" ,ylim=yrange)

plot(t, s1V, type="l", lty=1, col="black", xaxt="n", xlab="", ylab="", ylim=xrange-mean(xrange));
plot(t, s1H, type="l", lty=1, col="black", xaxt="n", xlab="", ylab="", ylim=yrange-mean(yrange));

plot(t, varicella*1000-s1V-trendV, type="l", lty=1, col="black", xlab="Time (Year)", ylab="", xaxp=c(2000,2010,10), ylim=xrange-mean(xrange));
axis(1,xaxp=c(0,10,10))

plot(t, HZ*10000-s1H-trendH, type="l", lty=1, col="black", xlab="Time (Year)", ylab="", xaxp=c(2000,2010,10),ylim=yrange-mean(yrange));
axis(1,xaxp=c(0,10,10))




