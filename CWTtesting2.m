close all
fs = 10e3;
t = 0:1/fs:2;
x = sin(2*pi*0.1*fs*t)+sin(2*pi*0.4*fs*t);
figure(1)
plot(x);
figure(2)
stft(x,fs,'Window',kaiser(256,5),'OverlapLength',220,'FFTLength',512);
