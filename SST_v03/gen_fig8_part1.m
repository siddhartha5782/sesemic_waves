%
% run SST on HZ. save the result for display
%
clear; close all;


varicella = 0;

    %% to test the signal with one or two components
sr = 52;


load 500.HZ.20120127.1mCohort.txt
x = 1000*X500_HZ_20120127_1mCohort(:,4)./X500_HZ_20120127_1mCohort(:,2);


    %% the fontsize of the final reports
fz = 20;


    %% time
t = [1/sr:1/sr:10]+2000;



[sst] = setup_param;
sst.symmetry = 0;
sst.freqrange.low = 0;
sst.freqrange.high = 26;
sst.TFR.FWHM = 0.3;
sst.TFR.alpha = 0.02;
sst.TFR.alpha = 1/25;
sst.wavelet.linear = 1;  
sst.reconstruction = 1;
sst.display.reconstruction = 0;
sst.display.AM = 0;
sst.display.iff = 0;
sst.extractcurv.no = 1;
sst.extractcurv.iff1 = 1.4;
sst.extractcurv.range = 10;
sst.extractcurv.lambda = 12;



[f, nup, n1, n2] = padsignal(x','symmetric');


sst.x = f;
sst.t = [1:length(sst.x)]/sr - n1/sr + 2000;

[rslt] = synchrosqueezing(sst);

x1hat = zeros(size(f));
allhat = zeros(size(f));

RR = 4;
TT = 20;	

for oo = 1: length(f)
    allhat(oo) = rslt.alpha*real(sum(rslt.stfd(oo, TT:end)))*2;
    x1hat(oo) = rslt.alpha*real(sum(rslt.stfd(oo, rslt.c1(oo)-RR:rslt.c1(oo)+RR)))*2;
end

trendhat = f - allhat;

S = 5;
if 1
    trendhat0 = trendhat;
    for pp = S:length(x1hat)-S
	if pp<18 
            trendhat(pp) = median(trendhat0(1:pp+S-1));
	elseif pp>length(x1hat)-S
            trendhat(pp) = median(trendhat0(pp:end));
	else
	    trendhat(pp) = median(trendhat0(pp-S:pp+S));
	end
    end
end

x1hat = x1hat(n1+1:n1+length(x));
trendhat = trendhat(n1+1:n1+length(x));
noisehat = f(n1+1:n1+length(x)) - x1hat - trendhat;



sst.freqrange.high = 3;
sst.display.SSTcurv = 0;
[rslt] = synchrosqueezing(sst);

subplot(rslt.h2)
set(gca,'fontsize',fz);
ylabel('');
xlabel('');
axis([1/sr+2000 length(x)/sr+2000 -inf inf]);


subplot(rslt.h1)
set(gca,'fontsize',fz);
plot(sst.t, sst.x, 'k');
axis([1/sr+2000 length(x)/sr+2000 -inf inf]);
I=legend(['$X_1$']);
set(I,'Interpreter','Latex');


figure;
subplot(4,1,1);
plot(t, x,'k','linewidth',1.2);
hold on;
plot(t, trendhat+x1hat,'r');
plot(t, trendhat,'b');
set(gca,'fontsize',fz); axis tight;
I = legend(['$Y$']);
set(I,'Interpreter','Latex');

subplot(4,1,2); 
plot(t, x1hat,'k','linewidth',1.2);
set(gca,'fontsize',fz); axis tight;
I = legend(['$\hat{s}$']);
set(I,'Interpreter','Latex');
%set(gca,'xtick',[]);


subplot(4,1,3);
plot(t, noisehat,'k');
axis tight; set(gca,'fontsize',fz);
I = legend(['$\hat{\Phi}$']);
set(I,'Interpreter','Latex');
%set(gca,'xtick',[]);

subplot(4,1,4);
plot(t, trendhat,'k','linewidth',1.2);
axis tight; set(gca,'fontsize',fz);
xlabel('time');
I = legend(['$\hat{T}$']);
set(I,'Interpreter','Latex');

stfd = rslt.stfd(n1+1:n1+520, :);
freq = rslt.freq';

save HZ_rslt stfd trendhat noisehat x1hat freq t x
