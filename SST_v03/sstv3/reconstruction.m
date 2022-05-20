function [recon] = reconstruction(sst, rslt, D)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
range = sst.extractcurv.range;
iff = sst.extractcurv.iff1;
freq = rslt.freq;
alpha = rslt.alpha;

if strcmp(sst.TFRtype,'CWT') & ~sst.TFR.linear
    Z1 = -floor(log2(range*iff)) + floor(sst.TFR.nvoice);
    Z2 = -floor(log2(range*iff)) + floor(sst.TFR.nvoice);
else
    Z1 = ceil(0.4/alpha);	%% 0.1 Hz up
    Z2 = ceil(0.4/alpha);	%% 0.1 Hz down
end

tmpc1 = rslt.c1-Z1;
tmpc2 = rslt.c1+Z2;
if min(tmpc1) < 1; Z1 = Z1 + min(tmpc1) - 1; end
if max(tmpc2) > length(freq); Z2 = Z2 - (max(tmpc2) - length(freq)); end


recon = zeros(1,size(rslt.stfd,1));

for jj = 1:D

    recon0 = zeros(1,size(rslt.stfd,1));

    if strcmp(sst.TFRtype,'CWT') & ~sst.TFR.linear

        deltaxi = log2(alpha);
	    
        for kk = 1:size(rslt.stfd,1)
            stfd2 = log(2)*deltaxi*stfd(kk, jj*rslt.c1(kk)-Z1:jj*rslt.c1(kk)+Z2) * diag( freq(jj*rslt.c1(kk)-Z1:jj*rslt.c1(kk)+Z2) );
       	    recon0(kk) = 2*sum( stfd2 );
        end

    else

        for kk = 1:size(rslt.stfd,1)
       	    stfd2 = alpha*rslt.stfd(kk, jj*rslt.c1(kk)-Z1: jj*rslt.c1(kk)+Z2);
            recon0(kk) = 2*sum( stfd2 );
        end

    end

    recon = recon + recon0;
end

