%
% run SST on T_{1,2,1} and T_{1,3,1} and save the results for display
%


fprintf('================================================\n');
fprintf('\tGenerate Table 2\n\n');


clear; close all;
initstate(2999);

denoise = 1;
sr = 100;
Tlen = 10;
    %% time
t = [1/sr:1/sr:Tlen]; t = t';


    %% number of realization
NN = 51;

CC = 0;

for ss0 = 1:2; for kkk = 2:3; for jj0 = 1:2;
CC = CC + 1;

    %% kkk=2 is X_2. kkk=3 is X_3

    %% trend type, 1 and 2
Ttype = jj0 ;

    %% noise level, 0.5 and 1
sigma0 = ss0/2 ; 

fprintf([num2str(CC),'/8   ']);
fprintf(['Run X',num2str(kkk),' with T',num2str(jj0),' sigma0 = ',num2str(sigma0),'\n']);



    %% load the noise generated in R for the fair comparison
load t4_sim1_noise.mat ;
load garch_noise.mat

    %% this is the seasonality   
    %% this is the first component
AM0 = (1+0.1*cos(t)) .* transpose(atan(fliplr(linspace(-10,13,length(t)))))/2+2;

    %% change it to get more simulations (not shown in the paper)
change = 1; 
AM1 = AM0.*change;
phase1 = t+0.1*sin(t);
x1 = AM1.*cos(2*pi*phase1);

    %% this is the second component
AM2 = 3.5*ones(size(t)); AM2(751:1000) = 2;
phase2 = 3.4*t-0.02*t.^(2.3);
x2 = AM2.*cos(2*pi*phase2);



    %% this is the trend
if Ttype == 1
    trend = 8*(1./(1+(t/5).^2)+exp(-t/10));
elseif Ttype == 2
    trend = 2*t+10*exp(-((t-4).^2)/6);
end  

    %% this is the clean signal
x = x1 + x2 + trend;


RR1 = 12;
RR2 = 10;
TT = 20;    



[sst] = setup_param;
sst.symmetry = 0;
sst.freqrange.low = 0;
sst.freqrange.high = 40;
sst.TFR.FWHM = 0.3;	
sst.TFR.alpha = 0.02;
sst.wavelet.linear = 1;  
sst.reconstruction = 1;
sst.display.reconstruction = 0;
sst.display.SSTcurv = 0;
sst.display.reconband = 1;
sst.display.RR1 = RR1;
sst.display.RR2 = RR2;
sst.display.AM = 0;
sst.display.iff = 0;
sst.extractcurv.no = 2;
sst.extractcurv.iff1 = 1.2;
sst.extractcurv.iff2 = 3;
sst.extractcurv.range = 30;
sst.extractcurv.lambda = 15;



runtime = zeros(NN,1);

    %% accuracy
err_s1 = zeros(NN,1);
err_s2 = zeros(NN,1);
err_s = zeros(NN,1);
err_arma = zeros(NN,1);
err_T = zeros(NN,1);

save_s1hat = zeros(Tlen*sr, NN);
save_s2hat = zeros(Tlen*sr, NN); 
save_noisehat = zeros(Tlen*sr, NN);
save_That = zeros(Tlen*sr, NN); 
save_noise = zeros(Tlen*sr, NN);

for ll = 1:NN

    initstate(2999+ll);

    if kkk == 2
        noise = sigma0*Rdata(:,ll+1);
        noise(1:500) = noise(1:500)*0.8;
    elseif kkk == 3
	noise = garch(:, ll);
        noise = 2*sigma0*noise;
    end

    x = x(1:Tlen*sr) ;
    noise = noise(1:Tlen*sr) ;

    [f, nup, n1, n2] = padsignal(x + noise,'symmetric');
    save_noise(:,ll) = noise;

    f = f(:);
    sst.x = f;
    sst.t = [1:length(sst.x)]/sr - n1/sr;

    t1 = tic;
    [rslt] = synchrosqueezing(sst);

    x1hatC = zeros(size(f));
    x2hatC = zeros(size(f));
    allhat = zeros(size(f));

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
        phi1hat = smooth(phi1hat0, 20, 'loess');
        phi2hat = smooth(phi2hat0, 20, 'loess');
        a1hat = smooth(a1hat0, 20, 'loess');
        a2hat = smooth(a2hat0, 20, 'loess');
	x1hat = a1hat .* cos(phi1hat);
	x2hat = a2hat .* cos(phi2hat);
    end

    x1hat = x1hat(n1+1:n1+length(x1));
    x2hat = x2hat(n1+1:n1+length(x1));
    trendhat = trendhat(n1+1:n1+length(x1));

    save_s1hat(:,ll) = x1hat;
    save_s2hat(:,ll) = x2hat;
    save_noisehat(:,ll) = noisehat;
    save_That(:,ll) = trendhat;

    runtime(ll) = toc(t1);
    err_s1(ll) = norm(x1-x1hat)./norm(x1);
    err_s2(ll) = norm(x2-x2hat)./norm(x2);
    err_s(ll) = norm(x1+x2-x1hat-x2hat)./norm(x1+x2);
    err_arma(ll) = norm(noise-noisehat)./norm(noise);
    err_T(ll) = norm(trend-trendhat)./norm(trend);

    if err_s1(ll) > 0.25
        subplot(rslt.h1); hold off;
        plot(t, x1); hold on; plot(t, x1hat, 'r'); axis tight;
        title(['err_s1 = ', num2str(err_s1(ll))]);
	subplot(rslt.h2); axis([0 10 1751 2000]);	
	%pause(0.1); savefig(['table2_',num2str(kkk),'_',num2str(ss0),'_',num2str(ll),'.png'],'png');
    end

    close all

end

fprintf(['\nerror (s1, s2, T, arma):\n',num2str(mean(err_s1)),'\\pm',num2str(std(err_s1)),'\n',num2str(mean(err_s2)),'\\pm',num2str(std(err_s2)),'\n', num2str(mean(err_T)),'\\pm',num2str(std(err_T)),'\n',num2str(mean(err_arma)),'\\pm',num2str(std(err_arma)),'\n\n']);
fprintf(['time :',num2str(mean(runtime)),'\\pm',num2str(std(runtime)),'\n\n']);

ll = find(err_arma == median(err_arma));

x1hat = save_s1hat(:,ll);
x2hat = save_s2hat(:,ll);
trendhat = save_That(:,ll);
noisehat = save_noisehat(:,ll);
noise = save_noise(:,ll);

    %% save the results so that we can generate the pretty pictures in R
stfd = rslt.stfd(n1+1:4:n1+Tlen*sr, :);
freq = rslt.freq';
xt = x + noise;

eval(['save rslt/table2_X',num2str(kkk),'_T',num2str(jj0),'_',num2str(ss0),'_rslt stfd trend trendhat noise noisehat x1 x1hat x2 x2hat freq t xt']);

eval(['save rslt/table2_X',num2str(kkk),'_T',num2str(jj0),'_',num2str(ss0),'_allsave save_s1hat save_s2hat save_That save_noisehat']);

end; end; end

