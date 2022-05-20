function [sst,f]  = Wavelet_sst(x,fs)
N = numel(x);
x = x(:)';
dt = 1/fs;
wavparam = 6;
thr = 1e-8;
nv = 32;
noct = floor(log2(N))-1;

% Create scale vector
na = noct*nv;
a0 = 2^(1/nv);
scales = a0.^(1:na);
N_Sc = numel(scales);

%Create frequency vector for CWT computation
w = (1:fix(N/2));
w = w.*((2.*pi)/N);
w = [0, w, -w(fix((N-1)/2):-1:1)];

% Compute FFT of the time series
xdft = fft(x);
[wift,dwift]  = sstwaveft(w,scales,wavparam);

%Obtain CWT coefficients and derivative
cwtcfs = ifft(repmat(xdft,N_Sc,1).*wift,[],2);
dcwtcfs = ifft(repmat(xdft,N_Sc,1).*dwift,[],2);


%Compute the phase transform
phasetf = imag(dcwtcfs./cwtcfs)./(2*pi);


% Threshold for synchrosqueezing
phasetf(abs(phasetf)<thr) = NaN;

% Create frequency vector for output
log2Nyquist = log2(1/(2*dt));
log2Fund = log2(1/(N*dt));
freq = 2.^linspace(log2Fund,log2Nyquist,na);
Tx = 1/nv*sstalgo(cwtcfs,phasetf,thr);
sst = Tx;
f = freq;
function [wft,dwft] = sstwaveft(w,scales,wavparam)
N_Sc = numel(scales);
N_w = numel(w);
wft = zeros(N_Sc,N_w);
cf = wavparam;
%amor wavelet        
for j = 1:N_Sc
    expnt = -(scales(j).*w - cf).^2/2.*(w > 0);
    wft(j,:) = exp(expnt).*(w > 0);
end

%derivative of CWT
wMatrix = repmat(w,N_Sc,1);
dwft = 1j*wMatrix.*wft;
function Tx = sstalgo(cwtcfs,phasetf,gamma)
M = size(cwtcfs,1);
N = size(cwtcfs,2);
log2Fund = log2(1/N);
log2Nyquist = log2(1/2);
iRow = real(1 + floor(M/(log2Nyquist-log2Fund)*(log2(phasetf)-log2Fund)));
idxphasetf = find(iRow>0 & iRow<=M & ~isnan(iRow));
idxcwtcfs = find(abs(cwtcfs)>gamma);
idx = intersect(idxphasetf,idxcwtcfs);
iCol = repmat(1:N,M,1);
Tx = accumarray([iRow(idx) iCol(idx)],cwtcfs(idx),size(cwtcfs));