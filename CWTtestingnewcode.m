close all;
Fs=44100;
tclip=10e-3;
nos=Fs*tclip;
tpoints=0:10e-3:nos;
x=cos(2*pi*500*tpoints)+cos(2*pi*700*tpoints)+cos(2*pi*2700*tpoints)+cos(2*pi*1700*tpoints);
figure(1)
plot(tpoints,x);
figure(2)
cwt(x,'amor',seconds(1/Fs));
colormap(pink)
[cfs,perd]=cwt(x,'amor',seconds(1/Fs));
cfs=abs(cfs);
per=seconds(perd);
pcolor(tpoints,per,cfs);
set(gca,'yscale','log');
shading interp
figure(3)
pcolor(tpoints,per,cfs);
shading interp
figure(4)
[sst,f]=wsst(x,tpoints);
mesh(tpoints,f,abs(sst))
axis tight
view(2)
shading interp
figure(5)
[s,w,n] = fsst(x);

mesh(n,w/pi,abs(s))
shading interp
axis tight
view(2)



figure(6)

stft(x,Fs,'Window',kaiser(256,5),'OverlapLength',220,'FFTLength',512);