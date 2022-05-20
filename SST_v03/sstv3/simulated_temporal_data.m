function [data] = simulated_temporal_data(data_param);
%
% This code will generate some standard simulated/real data for comparison.
% 
% INPUT: 
%    data_param.datatype : the type of signal you want to generate
%    data_param.noise.type : the type of noise you want to generate. 
%    data_param.noise.snrdb : the noise level, defined in the standard way
%    data_param.noise.ns : non-stationary noise. ** Simply by modulation **
%    data_param.Hz : the sampling rate of the signal
%    data_param.T : the length of the signal
%    data_param.nonunif : if you want to try nonuniform sampling, then the parameters for nonuniform setup is stored here.
%
% the supported data type: 
%     1/2: one FM or AMFM component 
%     3/4: two FM or AMFM components, no crossover 
%     5/6: two FM or AMFM components, crossover 
%     7/8: multiple components with close but different frequencies
%     9/10/11: shape function examples in the shape paper
%     12: strong AM decay
%     13: reaction in between
%     14: BM IF
%     31: real data, respiratory signal
%     32: real data, ECG signal
%
% the supported noise type:
%     1: white noise
%     2: color noise
%     3: ARMA
%     4: ARCH	%% the code is from UCSD_GARCH toolbox
%     5: poisson
%
% By Hau-tieng Wu, 2011-09-26  
% By Hau-tieng Wu, 2012-03-03  
%

T = data_param.T;
Hz = data_param.Hz;

T = 2^nextpow2(T*Hz)/Hz;
data.t = [1/Hz : 1/Hz : T];	%% the time
data.freqrange.low = 1;		%% default
data.freqrange.high = 12;    	%% default


%==================================================================
    %% generate noise
CH = ones(size(data.t));
if data_param.noise.ns; CH = 1+0.1*cos(pi*data.t); end

e = randn(size(data.t));

switch data_param.noise.type
    case 1
    	noise = e;
    case 2
        % color noise; pending
    case 3
	a = [1 0.5]; b = [1 0.4 0.3 0.2];
	noise = filter(b,a,e); 
    case 4
	[noise, H] = garchsimulate(length(e),[1, 0.2, 0.2, 0.3],1,2);  
	noise = noise';
    case 5
	noise = random('poiss',1,size(data.t)); 
    case 6
	px = [-600:600]/40;
        p = 1./(10*sqrt(2*pi))*exp((-(px-4).^2)./6) + 1./(10*sqrt(2*pi))*exp((-(px+4).^2)./6);
        noise = randpdf(p,px,[1,length(data.x)]);
end

noise = CH.*noise; 



%==================================================================
    %% generate signal
switch data_param.datatype  
case 1

    %% 1 FM component example
    data.x = cos(2*pi*(5*data.t+0.9*cos(data.t)));

    data.extractcurv.no = 1;
    data.extractcurv.iff1 = 5;

    sigma = sqrt( var(data.x)*10.^( -data_param.noise.snrdb/10 ) );

    data.clean.component = zeros(data.extractcurv.no, length(data.t));
    data.clean.AM 	 = zeros(data.extractcurv.no, length(data.t));
    data.clean.phase	 = zeros(data.extractcurv.no, length(data.t));
    data.clean.IF	 = zeros(data.extractcurv.no, length(data.t));

    data.clean.component = cos(2*pi*(5*data.t+0.9*cos(data.t)));
    data.clean.AM 	 = ones(size(data.t));
    data.clean.phase 	 = 5*data.t+0.9*cos(data.t);
    data.clean.IF	 = 5-0.9*sin(data.t);
   
case 2

    %% 1 AM-FM component example
    data.x = (1+0.2*sin(data.t).^2).*cos(2*pi*(5*data.t+0.9*cos(data.t)));

    data.extractcurv.no = 1;
    data.extractcurv.iff1 = 5;

    sigma = sqrt( var(data.x)*10.^( -data_param.noise.snrdb/10 ) );

    data.clean.component = zeros(data.extractcurv.no, length(data.t));
    data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
    data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
    data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

    data.clean.component = cos(2*pi*(5*data.t+0.9*cos(data.t)));
    data.clean.AM        = 1+0.2*sin(data.t).^2;
    data.clean.phase     = 5*data.t+0.9*cos(data.t);
    data.clean.IF        = 5-0.9*sin(data.t);


case 3
    
    %% 2 FM components example, no crossover
    data.x = cos(2*pi*(3*data.t-0.5*(data.t.^2/T)))+cos(2*pi*(5*data.t+cos(data.t)));

    data.extractcurv.no = 2;
    data.extractcurv.iff1 = 3;
    data.extractcurv.iff2 = 5;

    sigma = var(data.x)*10.^( -data_param.noise.snrdb/10 );

    data.clean.component = zeros(data.extractcurv.no, length(data.t));
    data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
    data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
    data.clean.IF        = zeros(data.extractcurv.no, length(data.t));
    
    data.clean.component(1,:) = cos(2*pi*(3*data.t-0.5*(data.t.^2/T)));
    data.clean.AM(1,:)        = ones(size(data.t));
    data.clean.phase(1,:)     = 3*data.t-0.5*(data.t.^2/T);
    data.clean.IF(1,:)        = 3*ones(size(data.t))-data.t/T;
    data.clean.component(2,:) = cos(2*pi*(5*data.t+cos(data.t)));
    data.clean.AM(2,:)        = ones(size(data.t));
    data.clean.phase(2,:)     = 5*data.t+cos(data.t);
    data.clean.IF(2,:)        = 5-sin(data.t);

case 4

    %% 2 AM-FM components example, no crossover
    data.x = (1+0.2*cos(data.t+1)).*cos(2*pi*(3*data.t-0.5*(data.t.^2/T)))+sqrt(0.8+data.t/T).*cos(2*pi*(5*data.t+cos(data.t)));

    data.extractcurv.no = 2;
    data.extractcurv.iff1 = 3;
    data.extractcurv.iff2 = 5;

    sigma = var(data.x)*10.^( -data_param.noise.snrdb/10 );

    data.clean.component = zeros(data.extractcurv.no, length(data.t));
    data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
    data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
    data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

    data.clean.component(1,:) = cos(2*pi*(3*data.t-0.5*(data.t.^2/T)));
    data.clean.AM(1,:)        = 1+0.2*cos(data.t+1);
    data.clean.phase(1,:)     = 3*data.t-0.5*(data.t.^2/T);
    data.clean.IF(1,:)        = 3*ones(size(data.t))-data.t/T;
    data.clean.component(2,:) = cos(2*pi*(5*data.t+cos(data.t)));
    data.clean.AM(2,:)        = sqrt(0.8+data.t/T);
    data.clean.phase(2,:)     = 5*data.t+cos(data.t);
    data.clean.IF(2,:)        = 5-sin(data.t);

case 5
    
    %% 2 FM components example, crossover
    data.x = cos(2*pi*(2*data.t+(data.t).^(1.6)))+cos(2*pi*(data.t*8-(data.t).^(1.4)));

    data.extractcurv.no = 2;
    data.extractcurv.iff1 = 4;
    data.extractcurv.iff2 = 6;

    sigma = var(data.x)*10.^( -data_param.noise.snrdb/10 );

    data.clean.component = zeros(data.extractcurv.no, length(data.t));
    data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
    data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
    data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

    data.clean.component(1,:) = cos(2*pi*(2*data.t+(data.t).^(1.6)));
    data.clean.AM(1,:)        = ones(size(data.t));
    data.clean.phase(1,:)     = 2*data.t+(data.t).^(1.6);
    data.clean.IF(1,:)        = 2+1.6*(data.t).^(0.6);
    data.clean.component(2,:) = cos(2*pi*(data.t*8-(data.t).^(1.4)));
    data.clean.AM(2,:)        = ones(size(data.t));
    data.clean.phase(2,:)     = data.t*8-(data.t).^(1.4);
    data.clean.IF(2,:)        = 8-1.4*(data.t).^(0.4);

case 6

    %% 2 AM-FM components example, crossover
    data.x = (1+0.2*cos(data.t+1)).*cos(2*pi*(2*data.t+(data.t).^(1.6)))+sqrt(0.8+data.t/T).*cos(2*pi*(data.t*8-(data.t).^(1.4)));

    data.extractcurv.no = 2;
    data.extractcurv.iff1 = 4;
    data.extractcurv.iff2 = 6;

    sigma = var(data.x)*10.^( -data_param.noise.snrdb/10 );

    data.clean.component = zeros(data.extractcurv.no, length(data.t));
    data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
    data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
    data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

    data.clean.component(1,:) = cos(2*pi*(2*data.t+(data.t).^(1.6)));
    data.clean.AM(1,:)        = 1+0.2*cos(data.t+1);
    data.clean.phase(1,:)     = 2*data.t+(data.t).^(1.6);
    data.clean.IF(1,:)        = 2+1.6*(data.t).^(0.6);
    data.clean.component(2,:) = cos(2*pi*(data.t*8-(data.t).^(1.4)));
    data.clean.AM(2,:)        = sqrt(0.8+data.t/T);
    data.clean.phase(2,:)     = data.t*8-(data.t).^(1.4);
    data.clean.IF(2,:)        = 8-1.4*(data.t).^(0.4);

case 7

    %% two close frequency chirps
    data.x = (1+0.2*cos(data.t+1)).*cos(2*pi*(data.t+(data.t).^(1.5)))+sqrt(0.8+data.t/T).*cos(2*pi*(data.t*2+(data.t).^(1.5)));

    data.extractcurv.no = 2;
    data.extractcurv.iff1 = 3;
    data.extractcurv.iff2 = 5;

    sigma = var(data.x)*10.^( -data_param.noise.snrdb/10 );
    data.freqrange.low = 0.5;

    data.clean.component = zeros(data.extractcurv.no, length(data.t));
    data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
    data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
    data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

    data.clean.component(1,:) = (1+0.2*cos(data.t+1)).*cos(2*pi*(data.t+(data.t).^(1.5)));   
    data.clean.AM(1,:)        = 1+0.2*cos(data.t+1);
    data.clean.phase(1,:)     = data.t+(data.t).^(1.5);
    data.clean.IF(1,:)        = 1+1.5*(data.t).^(0.5);
    data.clean.component(2,:) = sqrt(0.8+data.t/T).*cos(2*pi*(data.t*2+(data.t).^(1.5)));   
    data.clean.AM(2,:)        = sqrt(0.8+data.t/T);
    data.clean.phase(2,:)     = data.t*2+(data.t).^(1.5);
    data.clean.IF(2,:)        = 2+1.5*(data.t).^(0.5);

case 8

    %% multiple harmonic components example
    data.x = cos(2*pi*(9*data.t))+cos(2*pi*(10*data.t))+cos(2*pi*(11*data.t))+cos(2*pi*(12*data.t));

    data.extractcurv.no = 4;
    data.extractcurv.iff1 = 9;
    data.extractcurv.iff2 = 10;
    data.extractcurv.iff3 = 11;
    data.extractcurv.iff4 = 12;

    sigma = var(data.x)*10.^( -data_param.noise.snrdb/10 );
    data.freqrange.low = 5;
    data.freqrange.high = 18;

    data.clean.component = zeros(data.extractcurv.no, length(data.t));
    data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
    data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
    data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

    data.clean.component(1,:) = cos(2*pi*(9*data.t));
    data.clean.AM(1,:)        = ones(size(data.t));
    data.clean.phase(1,:)     = 9*data.t;
    data.clean.IF(1,:)        = 9*ones(size(data.t));
    data.clean.component(2,:) = cos(2*pi*(10*data.t));
    data.clean.AM(2,:)        = ones(size(data.t));
    data.clean.phase(2,:)     = 10*data.t;
    data.clean.IF(2,:)        = 10*ones(size(data.t));
    data.clean.component(3,:) = cos(2*pi*(11*data.t));
    data.clean.AM(3,:)        = ones(size(data.t));
    data.clean.phase(3,:)     = 11*data.t;
    data.clean.IF(3,:)        = 11*ones(size(data.t));
    data.clean.component(4,:) = cos(2*pi*(12*data.t));
    data.clean.AM(4,:)        = ones(size(data.t));
    data.clean.phase(4,:)     = 12*data.t;
    data.clean.IF(4,:)        = 12*ones(size(data.t));

%    %% nonuniform examples. Figure 1 in the STFT paper
%    if nargin < 6
%        data_param.nonunif.NONUNIF = 8; 	%% the non-uniformality
%        data_param.nonunif.SS = 10;		%% the value used to control the uniform sampling rate. For example, if Hz=100, data_param.nonunif.SS=25, it means the sampling rate is Hz/data_param.nonunif.SS=4Hz. 
%    end
%
%    data.components = 1;
%    origx = (2+0.2*cos(data.t)).*cos(2*pi*(3*data.t+0.5*cos(data.t)));    
%    data.trueif(1).trueif = 3-0.5*sin(data.t);
%    data.trueam(1).trueam = 2+0.2*cos(data.t);
%
%    data.freqrange.low = 1;
%    data.origx = origx;
%    sigma = var(origx)*10.^( -data_param.noise.snrdb/10 );
%    origx = origx + sigma*randn(size(origx));
%    idx = data_param.nonunif.SS : data_param.nonunif.SS : length(data.t);
%    nonunifidx = round(rand(1,length(idx))*data_param.nonunif.NONUNIF); 
%    nonunifidx = idx + nonunifidx;
%    [nonunifidx] = sort(nonunifidx);
%    
%    data.x = zeros(size(origx));
%    data.nonunifidx = nonunifidx;
%    nonunifidx = nonunifidx(nonunifidx<T*Hz+1);
%    data.x(nonunifidx) = origx(nonunifidx);
%    for vv = 1:length(nonunifidx)-1
%        data.x(nonunifidx(vv)) = data.x(nonunifidx(vv))*(nonunifidx(vv+1)-nonunifidx(vv));
%    end
%
%elseif data_param.datatype == 8
%
%    %% nonuniform examples. Figure 5 in the STFT paper
%    if nargin < 6
%        data_param.nonunif.NONUNIF = 20;    %% the non-uniformality
%        data_param.nonunif.SS = 25;         %% the value used to control the uniform sampling rate. For example, if Hz=100, data_param.nonunif.SS=25, it means the sampling rate is Hz/data_param.nonunif.SS=4Hz. 
%    end
%
%    data.components = 1;
%    origx = cos(2*pi*(5*data.t+0.1*cos(data.t)));
%    data.trueif(1).trueif = 5-0.1*sin(data.t);
%    data.trueam(1).trueam = ones(size(data.t));
%
%    data.freqrange.low = 1;
%    data.origx = origx;
%    sigma = var(origx)*10.^( -data_param.noise.snrdb/10 );
%    origx = origx + sigma*randn(size(origx));
%    idx = data_param.nonunif.SS : data_param.nonunif.SS : length(data.t);
%    nonunifidx = round(rand(1,length(idx))*data_param.nonunif.NONUNIF);
%    nonunifidx = idx + nonunifidx;
%    [nonunifidx] = sort(nonunifidx);
%
%    data.nonunifidx = nonunifidx;
%    data.x = zeros(size(origx));
%    data.x(nonunifidx(nonunifidx<T*Hz+1)) = origx(nonunifidx(nonunifidx<T*Hz+1));
%    for vv = 1:length(nonunifidx(nonunifidx<T*Hz+1))-1
%        data.x(nonunifidx(vv)) = data.x(nonunifidx(vv))*(nonunifidx(vv+1)-nonunifidx(vv));
%    end    

case {9,10,11}

    %% shape function examples
    %% generate 3 different wave-shape functions in the shape paper
    t0 = [1/100:1/100:1];

    A0 = 1.2; %0.8;
    phi0 = t0;
    x0 = cos(A0*cos(2*pi*phi0))-sin(2*pi*phi0).*sin(A0*cos(2*pi*phi0))./cos(2*pi*phi0);
    x02 = (cos(2*pi*phi0+pi/3));
    x03 = zeros(size(x02));
    x03(find(x02>0)) = x02(find(x02>0)).^2;
    shape2 = (x0-1.4*x03).*cos(2*pi*(phi0));

    A0 = 1.2;
    phi0 = t0+0.3*cos(t0);
    shape3 = cos(2*pi*phi0).*cos(A0*cos(2*pi*phi0))-sin(2*pi*phi0).*sin(A0*cos(2*pi*phi0))./cos(2*pi*phi0);

    A0 = 0.8;
    shape1 = cos(2*pi*phi0).*cos(A0*cos(2*pi*phi0))-sin(2*pi*phi0).*sin(A0*cos(2*pi*phi0))./cos(2*pi*phi0);

    shape1 = shape1 - mean(shape1);
    shape1 = shape1/norm(shape1);
    shape1 = shape1/max(shape1);
    shape2 = shape2 - mean(shape2);
    shape2 = shape2/norm(shape2);
    shape2 = shape2/max(shape2);
    shape3 = shape3 - mean(shape3);
    shape3 = shape3/norm(shape3);
    shape3 = shape3/max(shape3);

    am1 = 1+0.1*sin(data.t.^(1.1));
    %phi1 = 1.5*data.t+0.1*cos(data.t+1);
    phi1 = 1.5*data.t+0.5*data.t.^(1.5);
    trueif1 = 1.5-0.2*sin(data.t+1);

    am2 = sqrt(1+0.2*cos(data.t));
    phi2 = 5*data.t+0.6*cos(data.t);
    trueif2 = 4.5*(1-0.2*sin(data.t));

    am3 = 1;
    phi3 = 7*data.t+cos(data.t).^2;
    trueif3 = 6.5-cos(data.t).*sin(data.t);

    if data_param.datatype == 9	%% one component

        data.x = am2.*shape2(round((mod(phi2,1)*99)+1));

        data.extractcurv.no = 1;
        data.extractcurv.iff1 = 4.5;

        data.clean.component = zeros(data.extractcurv.no, length(data.t));
        data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
        data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
        data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

    	data.clean.component(1,:) = am2.*shape3(round((mod(phi2,1)*99)+1));
    	data.clean.AM(1,:)        = am2;
    	data.clean.phase(1,:)     = phi2;
    	data.clean.IF(1,:)        = trueif2;
	data.clean.shape(1,:)	  = shape3;

    elseif data_param.datatype == 10	%% two components

        data.x = am1.*shape3(round((mod(phi1,1)*99)+1)) + am2.*shape2(round((mod(phi2,1)*99)+1));

        data.extractcurv.no = 2;
        data.extractcurv.iff1 = 1.5;
        data.extractcurv.iff2 = 4.5;

        data.clean.component = zeros(data.extractcurv.no, length(data.t));
        data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
        data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
        data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

        data.clean.component(1,:) = am1.*shape2(round((mod(phi1,1)*99)+1));
        data.clean.AM(1,:)        = am1;
        data.clean.phase(1,:)     = phi1;
        data.clean.IF(1,:)        = trueif1;
        data.clean.shape(1,:)     = shape2;
        data.clean.component(2,:) = am2.*shape3(round((mod(phi2,1)*99)+1));
        data.clean.AM(2,:)        = am2;
        data.clean.phase(2,:)     = phi2;
        data.clean.IF(2,:)        = trueif2;
        data.clean.shape(2,:)     = shape3;

    elseif data_param.datatype == 11	%% three components 

        data.x = am1.*shape1(round((mod(phi1,1)*99)+1)) + am2.*shape2(round((mod(phi2,1)*99)+1)) + am3.*shape3(round((mod(phi3,1)*99)+1));

        data.extractcurv.no = 3;
        data.extractcurv.iff1 = 1.5;
        data.extractcurv.iff2 = 4.5;
        data.extractcurv.iff3 = 6.5;

        data.clean.component = zeros(data.extractcurv.no, length(data.t));
        data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
        data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
        data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

        data.clean.component(1,:) = am1.*shape1(round((mod(phi1,1)*99)+1))
        data.clean.AM(1,:)        = am1;
        data.clean.phase(1,:)     = phi1;
        data.clean.IF(1,:)        = trueif1;
        data.clean.shape(1,:)     = shape1;
        data.clean.component(2,:) = am2.*shape2(round((mod(phi2,1)*99)+1))
        data.clean.AM(2,:)        = am2;
        data.clean.phase(2,:)     = phi2;
        data.clean.IF(2,:)        = trueif2;
        data.clean.shape(2,:)     = shape2;
        data.clean.component(3,:) = am3.*shape3(round((mod(phi3,1)*99)+1))
        data.clean.AM(3,:)        = am3;
        data.clean.phase(3,:)     = phi3;
        data.clean.IF(3,:)        = trueif3;
        data.clean.shape(3,:)     = shape3;

    end

    sigma = sqrt( var(data.x)*10.^( -data_param.noise.snrdb/10 ) );
    data.freqrange.high = 15;


case 14

    dd = clock;
    dd = dd(end);
    initstate(dd);

    %% 1 BM IF component example
    e = randn(size(data.t));
    iff = abs(cumsum(e));
    phi = cumsum(iff)*(data.t(2)-data.t(1));
    data.x = cos(2*pi*phi);

    data.extractcurv.no = 1;
    data.extractcurv.iff1 = mean(iff);

    sigma = sqrt( var(data.x)*10.^( -data_param.noise.snrdb/10 ) );
    data.freqrange.low = 0.1;        
    data.freqrange.high = 50;       

    data.clean.component = zeros(data.extractcurv.no, length(data.t));
    data.clean.AM        = zeros(data.extractcurv.no, length(data.t));
    data.clean.phase     = zeros(data.extractcurv.no, length(data.t));
    data.clean.IF        = zeros(data.extractcurv.no, length(data.t));

    data.clean.component = data.x;
    data.clean.AM        = ones(size(data.t));
    data.clean.phase     = phi;
    data.clean.IF        = iff;
   

case 31
    
    %% real world example. The respiratory signal
    load('resp.mat');
    data.t = t;
    data.x = x;
    Hz = 20;
    Tno = 2^(nextpow2(data.t(end)*Hz)-1);
    data.t = data.t(1:Tno);
    data.x = data.x(1:Tno);
    data.extractcurv.no = 1;
    data.extractcurv.iff1 = 0.3;
    data.freqrange.low = 0.1;
    data.freqrange.high = 1;
    noise = zeros(size(data.x));
    sigma = 0;

case 32

    %% real world example. The ECG signal
    load('ecg.mat');
    data.t = t(1:4:40000);
    data.x = zeros(size(data.t));

    for mm = 1: 10000
	data.x(mm) = max( x((mm-1)*4+1: mm*4) );
    end

    Hz = 250;
    Tno = 2^(nextpow2(data.t(end)*Hz)-1);
    data.t = data.t(1:Tno);
    data.x = data.x(1:Tno);
    data.extractcurv.no = 1;
    data.extractcurv.iff1 = 1;
    data.freqrange.low = 0.1;
    data.freqrange.high = 20;
    e = randn(1,4000);
    a = [1 0.5]; b = [1 0.4 0.3 0.2];
    noise = filter(b,a,e);
    data.x(1001:5000) = data.x(1001:5000) + noise;
    noise = zeros(size(data.x));
    sigma = 0;

end


data.x = data.x + sigma*noise;
data.clean.noise = sigma*noise;
data.sigma = sigma;
