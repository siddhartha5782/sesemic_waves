#
# run gen_fig4_part1.R before runing this code
# display the TBATS results on Y_{1,2,1} and Y_{1,3,1} 
#


par(mfrow = c(4,2));
par(mar=c(0.1,2,1,0.8)+0.1, oma=c(2,1,1,1), cex=2)

########################### 

time = (1:length(exportdata2$x))/100;

    # the observation
plot(time, exportdata2$x, type="l", xaxt="n", xlab="", ylab="" );
plot(time, exportdata3$x, type="l", xaxt="n", xlab="", ylab="" );


    # \hat{x}_1 and x_1
yrange = range(exportdata2$s1); 
#y2range = range(exportdata2$s1hat);
#yrange = range(exportdata2$s1); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
plot(time, exportdata2$s1, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);
par(new=TRUE); plot(time, exportdata2$x1hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)

yrange = range(exportdata3$s1); 
#y2range = range(exportdata3$s1hat);
#yrange = range(exportdata3$s1); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
plot(time, exportdata3$s1, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);
par(new=TRUE); plot(time, exportdata3$x1hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)


    # \hat{x}_2 and x_2
yrange = range(exportdata2$s2); 
#y2range = range(exportdata2$s2hat);
#yrange = range(exportdata2$s2); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
plot(time, exportdata2$s2, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);
par(new=TRUE); plot(time, exportdata2$x2hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)

yrange = range(exportdata3$s2); 
#y2range = range(exportdata3$s2hat);
#yrange = range(exportdata3$s2); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
plot(time, exportdata3$s2, type="l", col="red", xaxt="n", xlab="", ylab="", ylim=yrange);
par(new=TRUE); plot(time, exportdata3$x2hat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)


    # \hat{T} and T
yrange = range(exportdata2$trend); 
#y2range = range(exportdata2$trendhat);
#yrange = range(exportdata2$trend); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
plot(time, exportdata2$trend, type="l", col="red", ylim=yrange);
par(new=TRUE); plot(time, exportdata2$trendhat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)

yrange = range(exportdata3$trend); 
#y2range = range(exportdata3$trendhat);
#yrange = range(exportdata3$trend); yrange[1] = max(c(y1range[1],y2range[1])); yrange[2] = max(c(y1range[2],y2range[2]));
plot(time, exportdata3$trend, type="l", col="red", ylim=yrange);
par(new=TRUE); plot(time, exportdata3$trendhat, type="l", col="black", xaxt="n", yaxt="n", xlab="", ylab="", ylim=yrange)

