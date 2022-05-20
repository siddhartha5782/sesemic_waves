function [stftcoeff, tfd_ytick] = STFourier(tmpdata)
%
% Modified Short time Fourier transform ver 0.2
% Modified from wavelab 85
%
%
% INPUT:
%   tmpdata.x        : signal
%   tmpdata.t        : time
%   tmpdata.TFR      : parameters used for STFT
% (OPTIONS)
%   tmpdata.parallel : use parallel computataion
%   tmpdata.debug    : debug message
%
% OUTPUT:
%   stftcoeff 	     : short time Fourier transform
%
% DEPENDENCY:
%
% by Hau-tieng Wu 2011-06-20 (hauwu@math.princeton.edu)
%


x  = tmpdata.x; 
t  = tmpdata.t - tmpdata.t(1);
dt = t(2)-t(1);
n  = length(x);

if tmpdata.debug; fprintf('(DEBUG) start STFT\n'); end

xhat = fft(tmpdata.x);
xhat( ceil(length(xhat)/2)+1:end ) = 0;
    %% frequency axis index
xi = ( [(0:(n/2)) (((-n/2)+1):-1)]./n )./dt;

    %% fix the stepping size in the frequency axis (for computational/visualization issue)
if 1/(2*t(end))<tmpdata.alpha; 
    step = tmpdata.alpha; 
else 
    if tmpdata.debug; fprintf('alpha is too small. Use the Nyquist rate\n\n'); end
    step = 1/(2*t(end));
end

    %% the frequency range in Fourier side we would compute in STFT
MAXFREQ = tmpdata.MAXFREQ;
MINFREQ = tmpdata.MINFREQ;

if MAXFREQ > 1./(2*dt)
    error('in STFT, can''t exceed the sampling rate');
end

    %% total discretization points in the frequency axis
nxi = floor(MAXFREQ/step) - floor(MINFREQ/step) + 1;

if tmpdata.debug; 
    fprintf(['(STFT) Only focus on the frequency band [',num2str(MAXFREQ),',',num2str(MINFREQ),' Hz.\n(STFT) Total ', num2str(nxi),' frequencies will be calculated\n\n']);
end


	%% calculate R\psi
xi0 = [-5:1/10000:5];
if strcmp(tmpdata.STFTwindow, 'gaussian');   %% Gaussian (not really wavelet)

    Rpsi = 1;	%% the window is fixed (for now)

elseif strcmp(tmpdata.STFTwindow, 'Cinfc');

	%% the window is $\exp\{ \frac{1}{|x|^2-1} \}$
    window = exp( 1./( (xi0./tmpdata.FWHM ).^2-1 ) );
    window( find( xi0 >= tmpdata.FWHM ) ) = 0;
    window( find( xi0 <= -tmpdata.FWHM ) ) = 0;
    Rpsi = sum(window)/10000;

end


    %% freq axis index
tfd_ytick = [1:1:nxi]*step + MINFREQ;
stftcoeff = zeros(n,nxi);


    %% start to evaluate STFT
for ii = 1 : nxi

	xi0 = tfd_ytick(ii);

	if strcmp(tmpdata.STFTwindow, 'gaussian')

	        %% this kernel has m_0=1
	    window = exp(-(xi-xi0).^2./tmpdata.FWHM)./(sqrt(pi*tmpdata.FWHM));

	elseif strcmp(tmpdata.STFTwindow, 'Cinfc');

            window = exp( 1./( ((xi-xi0)./tmpdata.FWHM).^2-1 ) )./Rpsi;
            window( find( xi >= xi0 + tmpdata.FWHM ) ) = 0;
            window( find( xi <= xi0 - tmpdata.FWHM ) ) = 0;
 %if ~mod(ii,50); plot(xi, window); title([num2str(ii),'/',num2str(nxi)]); pause; end
	end


        what = window .* xhat;
        w = ifft(what);
        stftcoeff(:,ii) = w';

end

if tmpdata.debug; fprintf('(DEBUG) STFT is done\n');

end
