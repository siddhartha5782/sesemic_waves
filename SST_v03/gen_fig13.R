# to test the signal X2 ** with outliers **
#
# save data in Matlab:
# save simu2_rslt stfd trend trendhat noise noisehat x1 x1hat x2 x2hat freq t x

# load data here


QQQ1 <-readMat("/Users/hauwu/Dropbox/tex20120723/_WORKING_seasonality/released/rslt/simuY2_outlier_rslt.mat")
#QQQ1 <-readMat("/Users/hau-tiengwu/Dropbox/tex20120723/_WORKING_seasonality/released/simuY2_outlier_rslt.mat")


par(mfrow = c(4,2));
layout(matrix(c(1,2,3,4,3,5,3,6), 4, 2, byrow = TRUE))
par(mar=c(0.1,2,1,0.8)+0.1, oma=c(2,1,1,1), cex=1.4)



########################### 
    # the observation
ss = dim(QQQ1$t); if(ss[1]<ss[2]) time=t(QQQ1$t) else time=QQQ1$t;
ss = dim(QQQ1$xt); if(ss[1]<ss[2]) X=t(QQQ1$xt) else X=QQQ1$xt;
plot(time, X, type="l", xaxt="n", xlab="", ylab="");




    # \hat{x}_1 and x_1
y1range = range(QQQ1$x1); y2range = range(QQQ1$x1hat);
yrange = range(QQQ1$x1); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ1$x1); if(ss[1]<ss[2]) x1=t(QQQ1$x1) else x1=QQQ1$x1;
#plot(time, x1, type="l", lty=2, col=rgb(r=0.5,g=0.5,b=0.5), xaxt="n", xlab="", ylab="", ylim=yrange);
plot(time, x1, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);

ss = dim(QQQ1$x1hat); if(ss[1]<ss[2]) x1hat=t(QQQ1$x1hat) else x1hat=QQQ1$x1hat;
par(new=TRUE); plot(time, x1hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)



    # the SST of X
stfdt = (1:250)/25;
image(stfdt,QQQ1$freq,abs(QQQ1$stfd),col=rev(gray((0:64)/64)),xaxt="n", xlab="", ylab="Hertz")
axis(1,xaxp=c(0,10,5))
axis(2,xaxp=c(0,5,5))



    # \hat{x}_2 and x_2
y1range = range(QQQ1$x2); y2range = range(QQQ1$x2hat);
yrange = range(QQQ1$x2); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ1$x2); if(ss[1]<ss[2]) x2=t(QQQ1$x2) else x2=QQQ1$x2;
#plot(time, x2, type="l", lty=2, col=rgb(r=0.5,g=0.5,b=0.5), xaxt="n", xlab="", ylab="", ylim=yranage
plot(time, x2, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);

ss = dim(QQQ1$x2hat); if(ss[1]<ss[2]) x2hat=t(QQQ1$x2hat) else x2hat=QQQ1$x2hat;
par(new=TRUE); plot(time, x2hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)


    # \hat{T} and T
y1range = range(QQQ1$trend); y2range = range(QQQ1$trendhat);
yrange = range(QQQ1$trend); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ1$trend); if(ss[1]<ss[2]) trend=t(QQQ1$trend) else trend=QQQ1$trend;
#plot(time, trend, type="l", lty=2, col=rgb(r=0.5,g=0.5,b=0.5), xaxt="n", ylim=yrange);
plot(time, trend, type="l", col="red", xaxt="n", ylim=yrange);

ss = dim(QQQ1$trendhat); if(ss[1]<ss[2]) trendhat=t(QQQ1$trendhat) else trendhat=QQQ1$trendhat;
par(new=TRUE); plot(time, trendhat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)



    # \hat{Y} and Y
y1range = range(QQQ1$noise); y2range = range(QQQ1$noisehat);
yrange = range(QQQ1$noise); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ1$noise); if(ss[1]<ss[2]) noise=t(QQQ1$noise) else noise=QQQ1$noise;
#yrange = range(noise+shift); #yrange[1] = -5;
#plot(time, noise+shift, type="l", lty="longdash", col=rgb(r=0.5,g=0.5,b=0.5),ylim=yrange);
plot(time, noise, type="l", col="red", ylim=yrange);

#ss = dim(QQQ1$noisehat); if(ss[1]<ss[2]) noisediff=t(QQQ1$noisehat)-t(QQQ1$noise) else noisediff=QQQ1$noisehat-QQQ1$noise;
ss = dim(QQQ1$noisehat); if(ss[1]<ss[2]) noisehat=t(QQQ1$noisehat) else noisehat=QQQ1$noisehat;
par(new=TRUE); plot(time, noisehat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="",ylim=yrange)

