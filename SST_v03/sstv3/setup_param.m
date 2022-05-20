function [sst] = sst_param
%
% This is aimed to setup the parameters for the SST. 
%
% 2012-03-03 Hau-tieng Wu.
%

%% default parameters

    initstate; 	
	%% initiate the random seed

    sst.cleanup = 0;

    sst.debug = 0;
        %% 0/1: debug mode

    sst.symmetry = 0;
        %% 0/1: mirror symmetrize the signal to reduce the boundary effect

    sst.TFRtype = 'CWT';	
	%% CWT/STFT: CWT means using SST based on CWT.

    sst.TFR.nvoice = 64; 	
    sst.TFR.oct = 1;
    sst.TFR.scale = 2;
	%% nvoice, oct, scale: these are the parameters for CWT

    sst.TFR.motherwavelet = 'Cinfc';	
	%% morlet/gaussian/meyer/BL3/Cinfc: Now only these 5 types of 
	%% mother wavelets are supported. (TODO more)
    sst.TFR.STFTwindow = 'Cinfc';	
	%% gaussian/Cinfc: only these 2 types are supported (TODO more)

    sst.TFR.debug = sst.debug;
    sst.TFR.CENTER = 1;		
    sst.TFR.FWHM = 0.3;
	%% CENTER, FWHM: The parameters of the Gaussian wavelet

    sst.TFR.quantile_gamma = 0.1;	
	%% use quantile valaue as the threshold (1e-8 is ad-hoc default value)

    sst.TFR.alpha = 0.05;	
	%% the frequency resolution 

    sst.TFR.measure = 0;		
	%% if you want to use "measure version" (for STFT version)

    sst.TFR.linear = 1;	
	%% 0/1: use log/linear-scaled freq axis in the TF representation
	%% only for CWT version

    sst.TFR.rescale = 1; 
	%% used for STFT-SST. To rescale the time to estimate the low frequency

    sst.extractcurv.no      = 1;	
	%% Set up the number of curves you want to extract.
    sst.extractcurv.iff1    = 2;	
	%% iff{i}: The visual guess of the IF of the i-th component
    sst.extractcurv.lambda  = 10;	
	%% The regularity penalty	
    sst.extractcurv.range   = 48;	
	%% The range of TFD you run the curve extraction
    
    sst.display.signal	    = 1;
    sst.display.TFD	    = 0;	%% 0/1: display time-freq domain
    sst.display.SST	    = 1;	%% 0/1: display the SST
    sst.display.SSTcurv	    = 1;	%% 0/1: display the SST superimposed by the extracted curves
    sst.display.component   = 1;	%% 0/1: display the reconstructed comp
    sst.display.AM 	    = 1;	%% 0/1: display the reconstructed AM
    sst.display.phase 	    = 1;	%% 0/1: display the reconstructed phase
    sst.display.iff 	    = 1;	%% 0/1: display the reconstructed IFF
    sst.display.dominant    = 0;	%% 0/1: display the "dominant freq"
    sst.display.mesh	    = 1;	%% 0/1: display the mesh plot of SST
    sst.display.resolution  = 1;	%% the higher the lower the resolution
    sst.display.color	    = 'gray';	%% the colormap	
    sst.display.log	    = 0;	%% log scale in the TFR/SST value
    sst.display.fontsize    = 10;
    sst.display.SSTenhance  = 0;	%% enhance the local information in SST
    sst.display.timescale   = 'sec';	%% sec/min/hr/OTHER: the time scale
    sst.display.quantile    = 0.995;	%% to prevent extra large value

    sst.reconstruction = sst.display.component;
	%% reconstruct the AM/components or not
    

%% frequently used parameters

    sst.TFR.linear 	    = 1;	
    sst.TFR.quantile_gamma  = 0.01;	
    sst.TFR.FWHM            = 0.3;
    sst.TFR.alpha 	    = 0.05;	
    sst.extractcurv.lambda  = 10;
    sst.extractcurv.range   = 128;
    sst.freqrange.low       = 0.5;
    sst.freqrange.high      = 10;
    sst.display.SSTcurv	    = 0;	
    sst.display.component   = 0;
    sst.display.AM  	    = 0;
    sst.display.reconband   = 0;
    sst.display.curvband    = 0;
    sst.display.phase       = 0;
    sst.display.iff         = 0;
    sst.display.mesh	    = 0;
    sst.display.log	    = 0;
    sst.reconstruction      = 0;

    sst.TFR.MAXFREQ 	    = 20;
    sst.TFR.MINFREQ 	    = sst.TFR.FWHM;
