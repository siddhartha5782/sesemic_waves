function RSxx = reassignSpectrum(Sxx, f, t, fcorr, tcorr, options)
%reassignSpectrum move spectral estimates to their center of energy
%
%   This function is for internal use only. It may be removed in the future.
%
%   Copyright 2016-2019 The MathWorks, Inc.
%#codegen

nf = numel(f);
nt = numel(t);

fmin = f(1);
fmax = f(end);

tmin = t(1);
tmax = t(end);

% compute the destination row for each spectral estimate
% allow cyclic frequency reassignment only if we have a full spectrum
if isscalar(options.nfft) && strcmp(options.range,'twosided')
    rowIdx = 1+mod(round((fcorr(:)-fmin)*(nf-1)/(fmax-fmin)),nf);
else
    rowIdx = 1+round((fcorr(:)-fmin)*(nf-1)/(fmax-fmin));
end

% compute the destination column for each spectral estimate
colIdx = 1+round((tcorr(:)-tmin)*(nt-1)/(tmax-tmin));

% reassign the estimates that fit inside the time-frequency matrix
Sxx = Sxx(:);
idx = find(rowIdx>=1 & rowIdx<=nf & colIdx>=1 & colIdx<=nt);
RSxx = accumarray([rowIdx(idx) colIdx(idx)], Sxx(idx), [nf nt]);
