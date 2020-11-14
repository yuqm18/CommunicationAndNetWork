close all;
clear;
L = 10000; % number of symbol

Rs = 1;   % symbol rate in MHz
Fc = 1e3;   % carrier frequency in MHz
Gc = 0.8;   % channel gain
fs = Fc*10; % sample rate of analog wave

% fr = [0:F0/L:F0/2-F0/2 flip(-F0/L:-F0/L:-F0/2)];
fr = [0:Fc/L:fs-Fc/L; -fs:Fc/L:-Fc/L]';
wr = 2*pi*fr;
tau = 1/Rs/2;
MF = sum(fs*tau*sinc(wr*tau/2/pi).*exp(1j*wr*tau/2).*(1+0.5*exp(1j*(wr+2*pi*Fc)*100/3*10^8)),2);


figure;
plot([real(MF) imag(MF)])
figure;
plot([real(ifft(MF)) imag(ifft(MF))])

data = rand(L,1)>0.25;
vol = data.*(3+1i)+(1-data).*(-1-1i);

