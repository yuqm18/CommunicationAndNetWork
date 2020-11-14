clear;close all
bits = 7;
HalfVolNum = 2^(bits-1)-1;
Vint = 1/(2^(bits-1)-0.5);

delta = 0.01*Vint;

Sigma0 = 10.^(-6:0.05:1);
for nn = 1:length(Sigma0)
    Sigma = Sigma0(nn);
    delta1 = min(delta,Sigma/10);
y = -20*Sigma:delta1:20*Sigma;
Vol = (-HalfVolNum-0.5:1:HalfVolNum+0.5)*Vint;

error = min((y-Vol').^2);
Pe = sum(error.*exp(-y.^2/2/Sigma^2)/sqrt(2*pi)/Sigma *delta1);
% plot(y,error)
SNR(nn) = 10*log10(Sigma^2/Pe);
end
plot(Sigma0,SNR)
hold on
plot(Sigma0(Sigma0<1),10*log10(Sigma0(Sigma0<1).^2*12/Vint^2))
hh = gca;
hh.XScale = 'log';
title('高斯分布信源均匀量化信噪比')
ylabel('SNR/dB')
xlabel('Sigma/V')
legend('SNR_{q}','\sigma ^2/(\Delta^2/12)','Location','SouthEast')