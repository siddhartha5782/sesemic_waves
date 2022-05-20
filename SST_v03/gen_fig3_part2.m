%
% Run SST on Y0 (two harmonic components)
% run gen_noise.R first in R
%

clear; close all;
initstate(2999);

    %% to test the signal with one or two components
denoise = 1;
simulation = 0;
sr = 100;

    %% number of realization
NN = 201; 

    %% the noise generated in R for the fair comparison
if simulation == 0
    load t4_sim0_noise.mat;
elseif simulation == 1
    load t4_sim1_noise.mat;
else
    error();
end

    %% the fontsize of the final reports
fz = 20;


    %% time
t = [1/sr:1/sr:10]; t = t';

    %% this is the first component
x1 = 2.5*cos(2*pi*t);

    %% this is the second component
x2 = 3*cos(2*pi*pi*t); 



    %% this is the trend
trend = 8*(1./(1+(t/5).^2)+exp(-t/10)); 


    %% this is the clean signal
x = x1 + x2 + trend;


[sst] = setup_param;
sst.symmetry = 0;
sst.freqrange.low = 0;
sst.freqrange.high = 40;
sst.TFR.FWHM = 0.2;
sst.TFR.alpha = 0.02;
sst.wavelet.linear = 1;  
sst.reconstruction = 1;
sst.display.reconstruction = 0;
sst.display.SSTcurv = 1;
sst.display.AM = 0;
sst.display.iff = 0;
sst.extractcurv.no = 2;
sst.extractcurv.iff1 = 1;
sst.extractcurv.iff2 = 3.1;
sst.extractcurv.range = 10;
sst.extractcurv.lambda = 20;



runtime = zeros(NN,1);

    %% accuracy
err_x1 = zeros(NN,1);
err_x2 = zeros(NN,1);
err_arma = zeros(NN,1);
err_T = zeros(NN,1);
err_all = zeros(NN,1);

save_s1hat = zeros(10*sr, NN);
save_s2hat = zeros(10*sr, NN); 
save_noisehat = zeros(10*sr, NN);
save_That = zeros(10*sr, NN); 
save_noise = zeros(10*sr, NN);


fprintf('%4d',1);
for ll = 1:NN

    initstate(2999+ll);

    noise = Rdata(:,ll+1);
    if simulation == 1
        noise(1:500) = noise(1:500)*0.8;
    end

    [f, nup, n1, n2] = padsignal(x + noise,'symmetric');

    f = f(:);
    sst.x = f;
    sst.t = [1:length(sst.x)]/sr - n1/sr;

    t1 = tic;
    [rslt] = synchrosqueezing(sst);

    x1hatC = zeros(size(f));
    x2hatC = zeros(size(f));
    allhat = zeros(size(f));

    RR1 = 10;
    RR2 = 10;
    TT = 10;	

    for oo = 1: length(f)
	allhat(oo) = rslt.alpha*real(sum(rslt.stfd(oo, TT:end)))*2;
	x1hatC(oo) = rslt.alpha*sum(rslt.stfd(oo, rslt.c1(oo)-RR1:rslt.c1(oo)+RR1))*2;
	x2hatC(oo) = rslt.alpha*sum(rslt.stfd(oo, rslt.c2(oo)-RR2:rslt.c2(oo)+RR2))*2;
    end

    trendhat = f - allhat;

    noisehat = f - real(x1hatC) - real(x2hatC) - trendhat;
    noisehat = noisehat(n1+1:n1+length(x1));

    if denoise
        trendhat0 = trendhat;
        phi1hat0 = phase(x1hatC);
        phi2hat0 = phase(x2hatC);
        a1hat0 = abs(x1hatC);
        a2hat0 = abs(x2hatC);
        for pp = 11:length(x1hatC)-11
            trendhat(pp) = median(trendhat0(pp-10:pp+10));
        end
        phi1hat = smooth(phi1hat0, 50, 'loess');
        phi2hat = smooth(phi2hat0, 50, 'loess');
        a1hat = smooth(a1hat0, 50, 'loess');
        a2hat = smooth(a2hat0, 50, 'loess');
        x1hat = a1hat .* cos(phi1hat);
        x2hat = a2hat .* cos(phi2hat);
    end


    x1hat = x1hat(n1+1:n1+length(x1));
    x2hat = x2hat(n1+1:n1+length(x1));
    trendhat = trendhat(n1+1:n1+length(x1));
    noisehat = f(n1+1:n1+length(x1)) - x1hat - x2hat - trendhat;

    save_s1hat(:,ll) = x1hat;
    save_s2hat(:,ll) = x2hat;
    save_noisehat(:,ll) = noisehat;
    save_That(:,ll) = trendhat;

    runtime(ll) = toc(t1);
    err_x1(ll) = norm(x1-x1hat)./norm(x1);
    err_x2(ll) = norm(x2-x2hat)./norm(x2);
    err_arma(ll) = norm(noise-noisehat)./norm(noise);
    err_T(ll) = norm(trend-trendhat)./norm(trend);
    err_all(ll) = err_x1(ll) + err_x2(ll) + err_arma(ll) + err_T(ll);

    subplot(rslt.h2); axis([0 10 1801 2000]);
    subplot(rslt.h1); axis([0 10 -inf inf]);

    fprintf('\b\b\b\b'); fprintf('%4d',ll);

end

fprintf(['\nerror (s1, s2, T, arma):\n',num2str(mean(err_x1)),'\\pm',num2str(std(err_x1)),'\n',num2str(mean(err_x2)),'\\pm',num2str(std(err_x2)),'\n', num2str(mean(err_T)),'\\pm',num2str(std(err_T)),'\n',num2str(mean(err_arma)),'\\pm',num2str(std(err_arma)),'\n\n']);
fprintf(['time :',num2str(mean(runtime)),'\\pm',num2str(std(runtime)),'\n\n']);

ll = find(err_arma == median(err_arma));
%ll = find(err_all == median(err_all));

x1hat = save_s1hat(:,ll);
x2hat = save_s2hat(:,ll);
trendhat = save_That(:,ll);
noisehat = save_noisehat(:,ll);
noise = Rdata(:,ll+1);


close all
    %% generate SST
[f, nup, n1, n2] = padsignal(x + noise, 'symmetric');
sst.x = f;
sst.t = [1:length(sst.x)]/sr - n1/sr;
sst.display.reconband = 1;
sst.display.RR1 = RR1;
sst.display.RR2 = RR2;
sst.display.SSTcurv = 1;

    %% this is for the display purpose
sst.freqrange.high = 5;
[rslt] = synchrosqueezing(sst);
close all;




    %% save the results so that we can draw the pretty pictures in R
stfd = rslt.stfd(n1+1:4:n1+1000, :);
freq = rslt.freq';
xt = x + noise;

if simulation == 0
    save simuY0_rslt stfd trend trendhat noise noisehat x1 x1hat x2 x2hat freq t xt
elseif simulation == 1
    save simuY1_rslt stfd trend trendhat noise noisehat x1 x1hat x2 x2hat freq t xt
end

