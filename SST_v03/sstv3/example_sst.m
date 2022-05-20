%
% This is the example for you to try the Synchrosqueezing code. 
% 
% The first part of this file is generating some simulated signal
% The second part of this file is setting the parameters used in the main
%	code synchrosqueezing.m. The meaning of each parameter can be found
%	in synchrosqueezing.m
%
% 2012-03-03 Hau-tieng Wu.
%

%clear;


data_param.datatype = 10; 
    %% the type of signal you want to generate. See below

data_param.noise.type = 3; 
    %% the type of noise you want to generate. See below

data_param.noise.snrdb = 0; 
    %% the noise level, defined in the standard way

data_param.noise.ns = 1;
    %% the noise is nonstationary. Simply by modulating the variance

data_param.Hz = 200; 
    %% the sampling rate of the signal

data_param.T = 10; 
    %% the length of the signal, in second

data_param.nonunif = 0; 
    %% if you want to try nonuniform sampling, then the parameters for nonuniform setup is stored here.

[data] = simulated_temporal_data(data_param);

%
% the supported data type:
%     1/2: one FM or AMFM component
%     3/4: two FM or AMFM components, no crossover
%     5/6: two FM or AMFM components, crossover
%     7/8: multiple components with close but different frequencies
%     9/10/11: shape function examples in the shape paper
%     12: strong AM decay
%     13: reaction in between
%     31: real data, respiratory signal
%     32: real data, ECG signal
%
% the supported noise type:
%     1: white noise
%     2: color noise
%     3: ARMA
%     4: ARCH
%     5: poisson


	%%
	%% PART II: Setup parameters for synchrosqueezing.m
	%%

[sst] = setup_param;
sst.t = data.t;
sst.x = data.x;
sst.TFRtype = 'CWT';
sst.freqrange = data.freqrange;
sst.extractcurv.no = data.extractcurv.no;
sst.freqrange.high = 15;
sst.freqrange.low = 0;
sst.TFR.alpha = 0.02;
sst.TFR.MAXFREQ = 60;
sst.TFR.MINFREQ = 0;
sst.TFR.FWHM = 0.3;

for ii = 1:sst.extractcurv.no
    eval(['sst.extractcurv.iff',num2str(ii),' = data.extractcurv.iff',num2str(ii),';']);
end
 
[rslt] = synchrosqueezing(sst);
