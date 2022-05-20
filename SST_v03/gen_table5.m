%
% generate table 5 (sensitivity issue of SST)
%

fprintf('================================================\n');
fprintf('\tGenerate Table 5\n\n');


clear; close all;
initstate(2999);

denoise = 1;
sr = 100;
Tlen = 9;
    %% time
t = [1/sr:1/sr:Tlen]; t = t';


    %% number of realization
NN = 201;


CC = 0;
for kk0 = 2:3; for jj0 = 1:2; for ss0 = 1:2

CC = CC + 1;
    %% to test the signal with two components
    %% simulation = 2 is Y_1; = 3 is Y_2
simulation = kk0 ;

    %% trend type, 1 and 2
Ttype = jj0 ;

    %% noise level, 0.5 and 1
sigma0 = ss0/2 ; 



fprintf([num2str(CC),'/8   ']);
fprintf(['Run sensitivity X',num2str(kk0),' with T',num2str(jj0),' sigma0 = ',num2str(sigma0),'\n']);

eval(['load rslt/table2_X',num2str(kk0),'_T',num2str(jj0),'_',num2str(ss0),'_allsave']);


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
AM2 = 3.5*ones(size(t)); AM2(751:sr*Tlen) = 2;
phase2 = 3.4*t-0.02*t.^(2.3);
x2 = AM2.*cos(2*pi*phase2);



    %% this is the trend
if Ttype == 1
    trend = 8*(1./(1+(t/5).^2)+exp(-t/10));
elseif Ttype == 2
    trend = 5*t+20*exp(-((t-4).^2)/3);
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
sst.TFR.FWHM = 0.2;	
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
sst.extractcurv.lambda = 20;



runtime = zeros(NN,1);

    %% accuracy
err2_s1 = zeros(NN,1);
err2_s2 = zeros(NN,1);
err2_s = zeros(NN,1);
err2_arma = zeros(NN,1);
err2_T = zeros(NN,1);

save2_s1hat = zeros(Tlen*sr, NN);
save2_s2hat = zeros(Tlen*sr, NN); 
save2_noisehat = zeros(Tlen*sr, NN);
save2_That = zeros(Tlen*sr, NN); 
save2_noise = zeros(Tlen*sr, NN);

for ll = 1:NN

    initstate(2999+ll);

    if simulation == 2
        noise = sigma0*Rdata(:,ll+1);
        noise(1:500) = noise(1:500)*0.8;
    elseif simulation == 3
	noise = garch(:, ll);
        noise = 2*sigma0*noise;
    end

    x = x(1:Tlen*sr) ;
    noise = noise(1:Tlen*sr) ;

    [f, nup, n1, n2] = padsignal(x + noise,'symmetric');

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
	x1hatC(oo) = rslt.alpha*(sum(rslt.stfd(oo, rslt.c1(oo)-RR1:rslt.c1(oo)+RR1)))*2;
	x2hatC(oo) = rslt.alpha*(sum(rslt.stfd(oo, rslt.c2(oo)-RR2:rslt.c2(oo)+RR2)))*2;
    end
 
    trendhat = f - allhat;
    noisehat = f - real(x1hatC) - real(x2hatC) - trendhat;
    noisehat = noisehat(n1+1:n1+length(x)) ;

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

    x1hat = x1hat(n1+1:n1+length(x));
    x2hat = x2hat(n1+1:n1+length(x));
    trendhat = trendhat(n1+1:n1+length(x));

    save2_s1hat(:,ll) = x1hat;
    save2_s2hat(:,ll) = x2hat;
    save2_noisehat(:,ll) = noisehat;
    save2_That(:,ll) = trendhat;
    save2_noise(:,ll) = noise;

    runtime(ll) = toc(t1);
    err2_s1(ll) = norm(save_s1hat(1:sr*Tlen,ll)-x1hat)./norm(x1hat);
    err2_s2(ll) = norm(save_s2hat(1:sr*Tlen,ll)-x2hat)./norm(x2hat);
    err2_s(ll) = norm(save_s1hat(1:sr*Tlen,ll)+save_s2hat(1:sr*Tlen,ll)-x1hat-x2hat)./norm(x1hat+x2hat);
    err2_arma(ll) = norm(save_noisehat(1:sr*Tlen,ll)-noisehat)./norm(noisehat);
    err2_T(ll) = norm(save_That(1:sr*Tlen,ll)-trendhat)./norm(trendhat);

    err_s1(ll) = norm(x1-x1hat)./norm(x1);
    err_s2(ll) = norm(x2-x2hat)./norm(x2);
    err_s(ll) = norm(x1+x2-x1hat-x2hat)./norm(x1+x2);
    err_arma(ll) = norm(noise-noisehat)./norm(noise);
    err_T(ll) = norm(trend-trendhat)./norm(trend);

    if err2_s1(ll) > 0.15
        subplot(rslt.h1); hold off; 
	plot(t, x1); hold on; plot(t, x1hat, 'r'); axis tight;
	plot(t, save_s1hat(1:sr*Tlen,ll), 'k'); axis tight;
        title(['err_s1 = ', num2str(err_s1(ll)),'; diff = ',num2str(err2_s1(ll))]);
        subplot(rslt.h2); axis([0 9 1751 2000]);
        %pause(0.1); savefig(['table5_',num2str(kk0),'_',num2str(ss0),'_',num2str(ll),'.png'],'png');
    end

    close all;
end


fprintf(['\nerror (s1, s2, T, arma):\n',num2str(mean(err_s1)),'\\pm',num2str(std(err_s1)),'\n',num2str(mean(err_s2)),'\\pm',num2str(std(err_s2)),'\n', num2str(mean(err_T)),'\\pm',num2str(std(err_T)),'\n',num2str(mean(err_arma)),'\\pm',num2str(std(err_arma)),'\n\n']);

fprintf(['\nestimation difference (s1, s2, T, arma):\n',num2str(mean(err2_s1)),'\\pm',num2str(std(err2_s1)),'\n',num2str(mean(err2_s2)),'\\pm',num2str(std(err2_s2)),'\n', num2str(mean(err2_T)),'\\pm',num2str(std(err2_T)),'\n',num2str(mean(err2_arma)),'\\pm',num2str(std(err2_arma)),'\n\n']);
fprintf(['time :',num2str(mean(runtime)),'\\pm',num2str(std(runtime)),'\n\n']);

ll = find(err_arma == median(err_arma));

x1hat = save2_s1hat(:, ll);
x2hat = save2_s2hat(:, ll);
trendhat = save2_That(:, ll);
noisehat = save2_noisehat(:, ll);
noise = save2_noise(:, ll);

x1hat10 = save_s1hat(1:sr*Tlen, ll);
x2hat10 = save_s2hat(1:sr*Tlen, ll);
trendhat10 = save_That(1:sr*Tlen, ll);
noisehat10 = save_noisehat(1:sr*Tlen, ll);

    %% save the results so that we can generate the pretty pictures in R
stfd = rslt.stfd(n1+1:4:n1+Tlen*sr, :);
freq = rslt.freq';
xt = x + noise;

eval(['save rslt/table5_X',num2str(kk0),'_T',num2str(jj0),'_',num2str(ss0),'_subI stfd trend trendhat trendhat10 noise noisehat noisehat10 x1 x1hat x1hat10 x2 x2hat x2hat10 freq t xt']);

end; end; end

