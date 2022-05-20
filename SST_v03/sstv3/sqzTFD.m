function [stfd, freq, alpha] = sqzTFD(squeezing)
%
% Synchro- and -squeezing the time-frequency domain ver 0.2
% called by synchrosqueezing.m
%
% INPUT:
%   squeezing.tftype	= CWT or STFT
%   squeezing.tfd       = time frequency domain, CWT or STFT
%   squeezing.w         = instantaneous freuquency information function
%   squeezing.TFR   = parameters used in continuous wavelet transform
%   squeezing.stft      = parameters used in short time Fourier transform
%   squeezing.freqrange = the frequency range displayed in the time frequency plane
% (OPTIONS)
%   squeezing.debug      = debug message
%
% OUTPUT:
%
% DEPENDENCY:
%   Wavelab (included), ContWavelet.m
%
% by Hau-tieng Wu 2011-06-20 (hauwu@math.princeton.edu)
% by Hau-tieng Wu 2012-03-03 (hauwu@math.princeton.edu)
%


if squeezing.debug; fprintf('(DEBUG:synchrosqueezing) working on SS...'); end

w = squeezing.w;

    %% plot absolute value of SST
w = abs(w);
tfd = squeezing.tfd;

    %% n is the number of sampline points; nscale is the number of scales in CWT or number of frequencies in STFT
[n, nscale] = size(tfd);


%=====================================================================
	%% CWT-based SST
if strcmp(squeezing.TFRtype, 'CWT')

	%% build up the frequency axis in SST
    if squeezing.TFR.linear

	    %% alpha is the resolution of the frequency axis
        alpha = squeezing.TFR.alpha;
        nalpha = floor((squeezing.freqrange.high-squeezing.freqrange.low)./alpha);
        stfd = zeros(n,nalpha);
        freq = ([1:1:nalpha])*alpha+squeezing.freqrange.low;
	nfreq = length(freq);
	
	if squeezing.debug;
    	    fprintf(['linear freq resolution is ',num2str(alpha),'\n']);   
	end

    else	%% log scale

	nfreq = nscale;
        stfd = zeros(size(tfd));
        freq = zeros(1, nfreq);
        freq(1) = squeezing.freqrange.low;

	    %% log scale resolution
        alpha = 2.^((log2(squeezing.freqrange.high/squeezing.freqrange.low))./(nscale-1));
        for k = 2: nscale; freq(k) = freq(k-1)*alpha; end

	if squeezing.debug;
    	    fprintf(['log2 freq resolution is ',num2str(alpha),'\n']);   
	end
	    %% (2^alpha)^{freq resolution} = 2;

    end
 

    for b = 1:n             %% Synchro-
        if squeezing.debug; if mod(b,floor(n/10))==0; fprintf('*'); end; end;

        for kscale = 1: nscale       %% -Squeezing

            qscale = squeezing.TFR.scale .* (2^(kscale/squeezing.TFR.nvoice));

		%% squeezing is really a reassignment according to the rule determined by w
	    if squeezing.TFR.linear	%% linear scale

                if (isfinite(w(b, kscale)) && (w(b, kscale)>0))
                    k = floor( ( w(b,kscale)-squeezing.freqrange.low )./ alpha )+1;

                    if (isfinite(k) && (k > 0) && (k < nfreq-1))
			ha = freq(k+1)-freq(k);
                    	stfd(b,k) = stfd(b,k) + log(2)*tfd(b,kscale)*sqrt(qscale)./ha/squeezing.TFR.nvoice;
                    end
                end

	    else  %% log scale

            	if (isfinite(w(b, kscale)) && (w(b, kscale)>0))
                    k = floor(( log2( w(b,kscale)/squeezing.freqrange.low ) )./ log2(alpha))+1;
                    if (isfinite(k) && (k > 0) && (k < nfreq))
			ha = freq(k+1)-freq(k);
                    	stfd(b,k) = stfd(b,k) + log(2)*tfd(b,kscale)*sqrt(qscale)./ha/squeezing.TFR.nvoice;
                    end
            	end

	    end
        end
    end
end   




%=====================================================================
        %% STFT-based SST
if strcmp(squeezing.TFRtype, 'STFT')

	%% setup the freq axis index in SST
    alpha = squeezing.TFR.alpha;
    nalpha = floor((squeezing.freqrange.high-squeezing.freqrange.low)./alpha);   
    stfd = zeros(n,nalpha);
    freq = ([1:1:nalpha])*alpha+squeezing.freqrange.low;

    for b = 1:n		%% Synchro-
        if squeezing.debug; if mod(b,floor(n/10))==0; fprintf('*'); end; end;
    
        for jj = 1:nscale		%% -Squeezing
            if isfinite(w(b,jj))
                k = floor((w(b,jj)-squeezing.freqrange.low)./alpha);
                %k = floor((+squeezing.tfd_ytic(jj) -w(b,jj)-squeezing.freqrange.low)./alpha);
                if (k>0) && (k<=nalpha);
		    if squeezing.TFR.measure
			    %% this is the measure version which sometimes help to increase the visualization. As is shown in SIAM paper
                        stfd(b,k) = stfd(b,k)+1;
		    else
			    %% this is the ordinary SST based on STFT
                        stfd(b,k) = stfd(b,k)+tfd(b,jj);
		    end
                end
            end
        end
    end
    
end

if squeezing.debug; fprintf('DONE!!\n'); end

