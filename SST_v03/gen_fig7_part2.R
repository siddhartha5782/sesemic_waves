#
# display SST on varicella
#
# save data in Matlab:
# save simu2_rslt stfd trend trendhat noise noisehat x1 x1hat x2 x2hat freq t x

# load data here


#QQQ1 <-readMat("/Users/hau-tiengwu/Dropbox/tex20120723/_WORKING_seasonality/fig/varicella_rslt.mat")
QQQ1 <-readMat("/Users/hauwu/Dropbox/tex20120723/_WORKING_seasonality/fig/varicella_rslt.mat")


par(mfrow = c(3,2));
layout(matrix(c(1,2,3,4,3,5), 3, 2, byrow = TRUE))
par(mar=c(0.1,2,1,0.8)+0.1, oma=c(2,1,1,1), cex=2)



########################### 
    # the observation
yrange = range(QQQ1$x);

ss = dim(QQQ1$t); if(ss[1]<ss[2]) time=t(QQQ1$t) else time=QQQ1$t;
ss = dim(QQQ1$x); if(ss[1]<ss[2]) X=t(QQQ1$x) else X=QQQ1$x;
plot(time, X, type="l", xaxt="n", xlab="", ylab="", ylim=yrange);
par(new=TRUE); plot(time, QQQ1$trendhat+QQQ1$x1hat, type="l", col="green", xaxt="n", yaxt="n", xlab="", ylab="" ,ylim=yrange)


    # \hat{T} and T
ss = dim(QQQ1$trendhat); if(ss[1]<ss[2]) trendhat=t(QQQ1$trendhat) else trendhat=QQQ1$trendhat;
plot(time, X, type="l", col="black", xaxt="n", xlab="", ylab="" )
par(new=TRUE); plot(time, trendhat, type="l", col="red", xaxt="n", yaxt="n", xlab="", ylab="" ,ylim=yrange)

    # the SST of X
image(time, QQQ1$freq, abs(QQQ1$stfd), col=rev(gray((0:64)/64)), xlab="", ylab="Hertz")
axis(1,xaxp=c(2000,2010,5))
axis(2,xaxp=c(0,5,5))

    # \hat{x}_1 and x_1
ss = dim(QQQ1$x1hat); if(ss[1]<ss[2]) x1hat=t(QQQ1$x1hat) else x1hat=QQQ1$x1hat;
plot(time, x1hat, type="l", col="black", xaxt="n", xlab="", ylab="", ylim=yrange-mean(yrange))


    # \hat{Y} and Y
ss = dim(QQQ1$noisehat); if(ss[1]<ss[2]) noise=t(QQQ1$noisehat) else noise=QQQ1$noisehat;
plot(time, noise, type="l",xaxp=c(2000,2010,10), ylim=yrange-mean(yrange));
axis(1,xaxp=c(0,10,10))

