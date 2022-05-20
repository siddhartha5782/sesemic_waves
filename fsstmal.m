function [sst,f,t] = fsstmal(x)

win = kaiser(min(256,length(x)),10);

Fs = 2*pi;

nwin = length(win);
nfft = nwin;
noverlap = nwin-1;

% Convert to column vectors
x = signal.internal.toColIfVect(x);
win = win(:);
Fs = Fs(1);

isSingle = isa(x,'single') || isa(win,'single');

if isSingle 
    convFact = 'single';
else
    convFact = 'double';
end


tout = cast(1:length(x),convFact);

if isSingle
    x = single(x);
    win1 = single(win);
    Fs = single(Fs);
else
    win1 = win;
    Fs = double(Fs);
end

% Pad the signal vector x
if mod(nwin,2)==1
    xp = [zeros((nwin-1)/2,1) ; x ; zeros((nwin-1)/2,1)];
else
    xp = [zeros((nwin)/2,1) ; x ; zeros((nwin-2)/2,1)];
end

nxp = length(xp);

% Place xp into columns for the STFT
xin = getSTFTColumns(xp,nxp,nwin,noverlap,Fs);

% Compute the STFT
[sstout,fout] = computeDFT(bsxfun(@times,win1,xin),nfft,Fs);
stftc = computeDFT(bsxfun(@times,dtwin(win1,Fs),xin),nfft,Fs);

% Compute the reassignment vector
fcorr = -imag(stftc./ sstout);
fcorr(~isfinite(fcorr)) = 0;
fcorr = bsxfun(@plus,fout,fcorr);
tcorr = bsxfun(@plus,tout,zeros(size(fcorr)));

% Multiply STFT by a linear phase shift to produce the modified STFT
m = floor(nwin/2);
inds = 0:nfft-1;
ez = exp(-1i*2*pi*m*inds/nfft)';
sstout = bsxfun(@times,sstout,ez);

% Reassign the modified STFT
options.range = 'twosided';
options.nfft = nfft;
sstout = reassignSpectrum(sstout, fout, tout, fcorr, tcorr, options);

% Reduce to one-sided spectra if the input is real, otherwise return a
% two-sided (centered) spectra.
if isreal(x)
    fout = psdfreqvec('npts',nfft,'Fs',Fs,'Range','half');
    sstout = sstout(1:length(fout),:);
else
    % Centered spectra
    sstout = centerest(sstout);
    fout = centerfreq(fout);
end
sst = sstout;
f = fout;
t = tout(:)';

     
