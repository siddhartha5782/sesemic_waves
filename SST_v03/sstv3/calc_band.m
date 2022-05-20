% called by synchrosqueezing.m and finalreport.m

    if strcmp(sst.TFRtype,'CWT') & ~sst.TFR.linear
        Z1 = -floor(log2(curv.range*curv.iff)) + floor(sst.TFR.nvoice);
        Z2 = -floor(log2(curv.range*curv.iff)) + floor(sst.TFR.nvoice);
    else
        Z1 = ceil(0.4/alpha);	%% 0.1 Hz up
        Z2 = ceil(0.4/alpha);	%% 0.1 Hz down
    end

    eval(['tmpc1 = rslt.c',num2str(qq),'-Z1;']);
    eval(['tmpc2 = rslt.c',num2str(qq),'+Z2;']);
    if min(tmpc1) < 1; Z1 = Z1 + min(tmpc1) - 1; end
    if max(tmpc2) > length(freq); Z2 = Z2 - (max(tmpc2) - length(freq)); end

