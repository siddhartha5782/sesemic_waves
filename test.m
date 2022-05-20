fs=120;
t=1/fs: 1/fs:4-1/fs;
s1=sin(2*pi*(10*t + 2*atan((2*t-2).^2)));
s2=sin(2*pi*(32*t+10*sin(t)));
s3=sin(2*pi*(44*t+10*sin(t)));
sig=s1+0.1*s2+0.01*s3;

[sst,f] = Wavelet_sst(sig,fs);
hp = pcolor(t,f,abs(sst));
hp.EdgeColor = 'none';
