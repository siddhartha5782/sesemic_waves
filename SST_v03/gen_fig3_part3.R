# display Figure 3 (TBATS and SST on two harmonic signals)
#
# save data in Matlab:
# save simu2_rslt stfd trend trendhat noise noisehat x1 x1hat x2 x2hat freq t x

# load data here


QQQ2 <-readMat("/Users/hauwu/Dropbox/tex20120723/_WORKING_seasonality/released/simuY0_rslt.mat")


par(mfrow = c(5,2));
par(mar=c(0.1,2,1,0.8)+0.1, oma=c(2,1,1,1), cex=1.2)



########################### 
    # the observation
time = t(QQQ1$t)
X = t(QQQ1$xt)
plot(time, X, type="l", xaxt="n", xlab="", ylab="", ylim=c(-8,28));

ss = dim(QQQ2$t); if(ss[1]<ss[2]) time=t(QQQ2$t) else time=QQQ2$t;
ss = dim(QQQ2$xt); if(ss[1]<ss[2]) X=t(QQQ2$xt) else X=QQQ2$xt;
plot(time, X, type="l", xaxt="n", xlab="", ylab="", ylim=c(-8,28));



    # \hat{x}_1 and x_1
y1range = range(QQQ1$x1); y2range = range(QQQ1$x1hat);
yrange = range(QQQ1$x1); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
x1 = t(QQQ1$x1)
plot(time, x1, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);

x1hat = t(QQQ1$x1hat)
par(new=TRUE); plot(time, x1hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)



y1range = range(QQQ2$x1); y2range = range(QQQ2$x1hat);
yrange = range(QQQ2$x1); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ2$x1); if(ss[1]<ss[2]) x1=t(QQQ2$x1) else x1=QQQ2$x1;
plot(time, x1, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);

ss = dim(QQQ2$x1hat); if(ss[1]<ss[2]) x1hat=t(QQQ2$x1hat) else x1hat=QQQ2$x1hat;
par(new=TRUE); plot(time, x1hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)


    # \hat{x}_2 and x_2
y1range = range(QQQ1$x2); y2range = range(QQQ1$x2hat);
yrange = range(QQQ1$x2); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
x2=t(QQQ1$x2)
plot(time, x2, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);
x2hat=t(QQQ1$x2hat)
par(new=TRUE); plot(time, x2hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)


y1range = range(QQQ2$x2); y2range = range(QQQ2$x2hat);
yrange = range(QQQ2$x2); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
x2 = t(QQQ2$x2)
plot(time, x2, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);

x2hat = t(QQQ2$x2hat)
par(new=TRUE); plot(time, x2hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)


    # \hat{T} and T
y1range = range(QQQ1$trend); y2range = range(QQQ1$trendhat);
yrange = range(QQQ1$trend); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
trend = t(QQQ1$trend)
plot(time, trend, type="l", col="red", xaxt="n", ylim=yrange);
trendhat = t(QQQ1$trendhat)
par(new=TRUE); plot(time, trendhat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)


y1range = range(QQQ2$trend); y2range = range(QQQ2$trendhat);
yrange = range(QQQ2$trend); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ2$trend); if(ss[1]<ss[2]) trend=t(QQQ2$trend) else trend=QQQ2$trend;
plot(time, trend, type="l", col="red", xaxt="n", ylim=yrange);

ss = dim(QQQ2$trendhat); if(ss[1]<ss[2]) trendhat=t(QQQ2$trendhat) else trendhat=QQQ2$trendhat;
par(new=TRUE); plot(time, trendhat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)


    # \hat{T} and T
noise = QQQ1$xt-QQQ1$x1-QQQ1$x2-QQQ1$trend;
y1range = range(noise); y2range = range(QQQ1$noisehat);
yrange = range(noise); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
noise = t(noise)
plot(time, noise, type="l", col="red", ylim=yrange);

noisehat = t(QQQ1$noisehat)
par(new=TRUE); plot(time, noisehat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange )


y1range = range(QQQ2$noise); y2range = range(QQQ2$noisehat);
yrange = range(QQQ2$noise); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ2$noise); if(ss[1]<ss[2]) noise=t(QQQ2$noise) else noise=QQQ2$noise;
plot(time, noise, type="l", col="red", ylim=yrange);

ss = dim(QQQ2$noisehat); if(ss[1]<ss[2]) noisehat = t(QQQ2$noisehat) else noisehat = QQQ2$noisehat;
par(new=TRUE); plot(time, noisehat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange )


