# for Figure 2. Display the SST result of the clean signal
#
# save data in Matlab:
# save simu2_rslt stfd trend trendhat noise noisehat x1 x1hat x2 x2hat freq t x

# load data here

QQQ1 <-readMat("/Users/hauwu/Desktop/code20111220/water/released/simu_clean_rslt.mat")


par(mfrow = c(3,2));
layout(matrix(c(1,2,1,3,1,4), 3, 2, byrow = TRUE))
par(mar=c(0.1,2,1,0.8)+0.1, oma=c(2,1,1,1), cex=2)


########################### 
    # the observation


    # the SST of X
image(t(QQQ1$stfdt),QQQ1$freq,abs(QQQ1$stfd),col=rev(gray((0:64)/64)),xaxt="n", xlab="", ylab="Hertz")
axis(1,xaxp=c(0,10,5))
axis(2,xaxp=c(0,5,5))


    # \hat{x}_1 and x_1
y1range = range(QQQ1$x1); y2range = range(QQQ1$x1hat);
yrange = range(QQQ1$x1); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ1$x1); if(ss[1]<ss[2]) x1=t(QQQ1$x1) else x1=QQQ1$x1;
plot(time, x1, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);

ss = dim(QQQ1$x1hat); if(ss[1]<ss[2]) x1hat=t(QQQ1$x1hat) else x1hat=QQQ1$x1hat;
par(new=TRUE); plot(time, x1hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)



    # \hat{x}_2 and x_2
y1range = range(QQQ1$x2); y2range = range(QQQ1$x2hat);
yrange = range(QQQ1$x2); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ1$x2); if(ss[1]<ss[2]) x2=t(QQQ1$x2) else x2=QQQ1$x2;
plot(time, x2, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);

ss = dim(QQQ1$x2hat); if(ss[1]<ss[2]) x2hat=t(QQQ1$x2hat) else x2hat=QQQ1$x2hat;
par(new=TRUE); plot(time, x2hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)


    # \hat{T} and T
y1range = range(QQQ1$trend); y2range = range(QQQ1$trendhat);
yrange = range(QQQ1$trend); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
ss = dim(QQQ1$trend); if(ss[1]<ss[2]) trend=t(QQQ1$trend) else trend=QQQ1$trend;
#plot(time, trend, type="l", lty=2, col=rgb(r=0.5,g=0.5,b=0.5), ylim=yrange, xlab="Time (Second)");
plot(time, trend, type="l", col="red", ylim=yrange, xlab="Time (Second)");

ss = dim(QQQ1$trendhat); if(ss[1]<ss[2]) trendhat=t(QQQ1$trendhat) else trendhat=QQQ1$trendhat;
par(new=TRUE); plot(time, trendhat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)

