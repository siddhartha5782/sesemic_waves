function [cwtcoeff, ftd_ytic, Rpsi] = ContWavelet(cwt)
%
% Continuous wavelet transform ver 0.1
% Modified from wavelab 85
%
% INPUT:
%   cwt.xhat     : signal
%   cwt.t        : time
%   cwt.oct	 : the initial octave
%   cwt.scale	 : the log scale
%   cwt.nvoice	 : the number of equally spaced bins in each log scale
%   cwt.motehrwavelet   : the mother wavelet used in CWT
%   cwt.minus_i_partial_b : evaluate -i\partial_bW_f(a,b)
%   cwt.CENTER	 : the center of the Gaussian function if the motherwavelet is Gaussian
%   cwt.FWHM  	 : the FWHM of the Gaussian fucntion if the motherwavelet is Gaussian
% (OPTIONS)
%   cwt.parallel : use parallel computataion
%   cwt.debug    : debug message
%
% OUTPUT:
%   cwtcoeff	 : Continuous wavelet transform on log scale
%
% DEPENDENCY:
%
% by Hau-tieng Wu 2011-06-20 (hauwu@math.princeton.edu)
%


xhat = cwt.xhat; 
t = cwt.t;

    %% sampling interval
dt = t(2)-t(1);
n = length(xhat);

if cwt.debug; fprintf('working on CWT\n'); end

    %% the frequency index of the fft
xi = [(0:(n/2)) (((-n/2)+1):-1)];

    %% how many octave used in the CWT
noctave = floor(log2(n)) - cwt.oct;

    %% nvoice = how many equally spaced interval inside an octave. In log2 scale
cwtcoeff = zeros(n, cwt.nvoice .* noctave);

    %% index of the scale
kscale = 1;

    %% the first scale in CWT
scale = cwt.scale;

    %% the index of the scale axis in CWT
ftd_ytic = zeros(1, cwt.nvoice .* noctave);

    %% build up the index of the scale axis in CWT
for jj = 1 : cwt.nvoice .* noctave
    ftd_ytic(jj) = scale .* (2^(jj/cwt.nvoice));
end


    %% perform the CWT here
for jo = 1:noctave	% # of scales
    for jv = 1:cwt.nvoice
        qscale = scale .* (2^(jv/cwt.nvoice));
	
        omega =  xi ./ qscale ;            

	    %% different mother wavelets with proper scale
	if strcmp(cwt.motherwavelet, 'morlet'); 	%% Morlet

    	    window = 4*sqrt(pi)*exp(-4*(omega-0.69*pi).^2)-4.89098d-4*4*sqrt(pi)*exp(-4*omega.^2);  

	elseif strcmp(cwt.motherwavelet, 'gaussian');	%% Gaussian (not really wavelet)

	    if ~isfield(cwt,'CENTER') | ~isfield(cwt,'FWHM')
		error('You should assign CENTER and FWHM if you want to use Gaussian function as the mother wavelet');
	    end

	    psihat = @(f) exp( -log(2)*( 2*(f-cwt.CENTER)./cwt.FWHM ).^2 );
	    %window = exp(-log(2)*(2*(omega-cwt.CENTER)./cwt.FWHM).^2);
	    window = psihat(omega);

	elseif strcmp(cwt.motherwavelet, 'Cinfc');	

	    window = exp( 1./( ((omega-cwt.CENTER)./cwt.FWHM).^2-1 ) );
	    window( find( omega >= (cwt.CENTER+cwt.FWHM) ) ) = 0;
	    window( find( omega <= (cwt.CENTER-cwt.FWHM) ) ) = 0;

	elseif strcmp(cwt.motherwavelet, 'meyer');	%% Meyer

            window = zeros(size(omega));
            int1 = find((omega>=5./8*0.69*pi)&(omega<0.69*pi));
            int2 = find((omega>=0.69*pi)&(omega<7./4*0.69*pi));
            window(int1) = sin(pi/2*meyeraux((omega(int1)-5./8*0.69*pi)/(3./8*0.69*pi)));
            window(int2) = cos(pi/2*meyeraux((omega(int2)-0.69*pi)/(3./4*0.69*pi)));

	elseif strcmp(cwt.motherwavelet, 'BL3');	%% B-L 3

            phihat = (2*pi)^(-0.5)*(sin(omega/4)./(omega/4)).^4; phihat(1) = (2*pi)^(-0.5);
            aux1 = 151./315 + 397./840*cos(omega/2) + 1./21*cos(omega) + 1./2520*cos(3*omega/2);
            phisharphat = phihat.*(aux1.^(-0.5));

            aux2 = 151./315 - 397./840*cos(omega/2) + 1./21*cos(omega) - 1./2520*cos(3*omega/2);
            aux3 = 151./315 + 397./840*cos(omega) + 1./21*cos(2*omega) + 1./2520*cos(3*omega);
            msharphat = sin(omega/4).^4.*(aux2.^(0.5)).*(aux3.^(-0.5));
            window = phisharphat.*msharphat.*exp(i*omega/2).*(omega>=0);

	else

	    error('You should choose a supported motherwave'); 

	end

        window = window ./ sqrt(qscale);

	if cwt.minus_i_partial_b
		%% fourier side implementation of partial_bW_f(a,b)
            what = (xi./(n*dt)) .* window .* xhat;
	else
		%% fourier side implementation of W_f(a,b)
            what = window .* xhat;
	end

        w = ifft(what);
        cwtcoeff(:,kscale) = transpose(w);
	    %% one nvoice
        kscale = kscale+1;

    end

	%% one octave
    scale = scale .* 2;

end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %% calculate the constant for reconstruction
    %% TODO: calculate Rpsi for other mother wavelets
xi = [0.05:1/10000:10];
if strcmp(cwt.motherwavelet, 'gaussian');   %% Gaussian (not really wavelet)

    window = exp(-log(2)*(2*(xi-cwt.CENTER)./cwt.FWHM).^2);
    Rpsi = sum(window./xi)/10000;

elseif strcmp(cwt.motherwavelet, 'Cinfc');

    window = exp( 1./( ((xi-cwt.CENTER)./cwt.FWHM).^2-1 ) );
    window( find( xi >= (cwt.CENTER+cwt.FWHM) ) ) = 0;
    window( find( xi <= (cwt.CENTER-cwt.FWHM) ) ) = 0;
    Rpsi = sum(window./xi)/10000;

end

    %% this normalization is necessary for reconstruction
cwtcoeff = cwtcoeff ./ Rpsi;


if cwt.debug; fprintf('CWT is done\n\n'); end
end

