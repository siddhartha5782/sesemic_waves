%
% for Figure 2. SST on the clean signal
%

clear; close all;
initstate(1999);


    %% to test the signal with one or two components
denoise = 1;
sr = 100;


    %% number of realization
NN = 1; 


    %% the fontsize of the final reports
fz = 28;



    %% time
t = [1/sr:1/sr:10];

    %% this is the seasonality   
    %% this is the first component
AM0 = (1+0.1*cos(t)) .* atan(fliplr(linspace(-10,13,length(t))))/2+2;
change = 1; 
AM1 = AM0.*change;
phase1 = t+0.1*sin(t)-1;
x1 = AM1.*cos(2*pi*phase1);

    %% this is the second component
AM2 = 3.5*ones(size(t)); AM2(751:1000) = 2;
phase2 = 3.4*t-0.02*t.^(2.3);
x2 = AM2.*cos(2*pi*phase2);



    %% this is the trend
trend = 8*(1./(1+(t/5).^2)+exp(-t/10));


    %% this is the clean signal
x = x1 + x2 + trend;


[sst] = setup_param;
sst.symmetry = 0;
sst.freqrange.low = 0;
sst.freqrange.high = 40;
sst.TFR.FWHM = 0.3;
sst.TFR.alpha = 0.02;
sst.wavelet.linear = 1;  
sst.reconstruction = 1;
sst.display.reconstruction = 0;
sst.display.AM = 0;
sst.display.iff = 0;
sst.extractcurv.no = 2;
sst.extractcurv.iff1 = 1;
sst.extractcurv.iff2 = 3;
sst.extractcurv.range = 35;
sst.extractcurv.lambda = 15;




    [f, nup, n1, n2] = padsignal(x,'symmetric');

    f = f(:);
    sst.x = f;
    sst.t = [1:length(sst.x)]/sr - n1/sr;

    t1 = tic;
    [rslt] = synchrosqueezing(sst);

    x1hat = zeros(size(x1'));
    x2hat = zeros(size(x1'));
    allhat = zeros(size(x1'));

    RR = 10;
    TT = 20;	%% 21 previously

    for oo = 1: length(x1)
	allhat(oo) = rslt.alpha*real(sum(rslt.stfd(n1+oo, TT:end)))*2;
	x1hat(oo) = rslt.alpha*real(sum(rslt.stfd(n1+oo, rslt.c1(n1+oo)-RR:rslt.c1(n1+oo)+RR)))*2;
	x2hat(oo) = rslt.alpha*real(sum(rslt.stfd(n1+oo, rslt.c2(n1+oo)-RR:rslt.c2(n1+oo)+RR)))*2;
    end

    trendhat = f(n1+1:n1+length(x1)) - allhat;

    if denoise
	trendhat0 = trendhat;
	for pp = 11:length(x1hat)-11
	    trendhat(pp) = median(trendhat0(pp-10:pp+10));
	end
    end


    %% for the display purpose
sst.freqrange.high = 5;
sst.display.SSTcurv = 0;
[rslt] = synchrosqueezing(sst);

close all;

    %% save the results for R 
stfd = rslt.stfd(n1+1:4:n1+1000, :);
freq = rslt.freq';
stfdt = t(1:4:end);

save simu_clean_rslt stfd stfdt trend trendhat x1 x1hat x2 x2hat freq t x
