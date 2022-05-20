load quadchirp;
n=1:1000;
x=sin(2*pi*n/5);
figure
[sst,f] = wsstmal(x);
plot(n,sst);
%figure()
%hp = pcolor(tquad,f,abs(sst));
%hp.EdgeColor = 'none';
%title('Wavelet Synchrosqueezed Transform');
%xlabel('Time'); ylabel('Hz');