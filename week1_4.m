T=3600;
f=(0:0.01:4)';
B=0.1;
fs = 20;
N0_0 = 0.1:0.05:1.6;%单边功率谱密度
rng(2);
ShowPSD = 0;

for repN0 = 1:length(N0_0)
    for repTrial = 1:10
        N0 = N0_0(repN0);

        x=rand(T,1);
        y=zeros(T,1);
        y(1)=1;
        for i=2:T
            if(y(i-1)==1)
                y(i)=(x(i)>.98)*2+(x(i)<=.98)*1;
            end
            if(y(i-1)==2)
                y(i)=(x(i)>.6)*-1+(x(i)<=.6)*2;
            end
            if(y(i-1)==0)
                y(i)=(x(i)>.6)*1+(x(i)<=.6)*0;
            end
            if(y(i-1)==-1)
                y(i)=(x(i)>.95)*0+(x(i)<=.95)*-1;
            end

        end
        y=y.*(y~=2);

        % y=y-mean(y);

        %y=sin(2*pi*0.2*(1:T)');%test


        t = (1/fs:1/fs:T)';
        y0 = y;
        y = y(ceil(t));

        % k=(min(f)-B/2)*T:(max(f)+B/2)*T;
        % % Xf=(1/T*(y.*exp(-2i*pi*t'.*k/T)))/fs;
        % % Xf=abs(sum(Xf)).^2;
        % 
        % Xf = abs(fft(y)/fs/T).^2;
        % Xf = Xf(mod(k,length(Xf))+1)';
        % 
        % k1=(k<=(f'+B/2)*T)&(k>=(f'-B/2)*T);
        % Xf1=k1.*Xf;
        % Xf1=sum(Xf1,2);

        if ShowPSD
        figure(2)
        Xf1 = PowerSpectralDensity(y,B,T,f,fs);

        plot(f,Xf1)
        xlabel('频率/Hz')
        ylabel('功率谱/Hz^{-1}')
        title('s_1')
        end
        figure(3)
        r = y+random('normal',0,N0/2,size(y));

        rx = filter(0.3,[1 -exp(-10/3/fs)],r);
        
        sys = zpk([],-1/0.3,1/0.3);
        rx = lsim(sys,r,t-t(1));

%         plot(rx)

        % 方法一，单独考虑每一个时刻
        rxs = rx(fs:fs:end);    %sampling
        x1 = (N0*log(4)-1)/2;   % -1 0 分界
        x2 = (-N0*log(10)+1)/2; % 0 1  分界
        x3 = -N0/2*log(5/2);    % -1 1 分界
        % x1==x2 时,x1 == x3
        if x1<x3
        rxs(rxs<x1) = -1;
        rxs(rxs<=x2&rxs>=x1) = 0;
        rxs(rxs>x2) = 1;
        else 
            rxs(rxs<x3) = -1;
            rxs(rxs>=x3) = 1;
        end

        err = abs(rxs-y0);
        ErrTrial(repTrial) = sum(err)/T;
        
        % 方法2 联合概率
        % dim 1: previous input, ax1
        % dim 2: current input, ax2
        % dim 3: next input, ax3
        Tm = [0.98 0  0.4 0;
            0 0.95 0 0.4;
            0 0.05 0.6 0;
            0.02 0 0 0.6]';
        T1 = reshape(Tm,4,1,4);
        T2 = Tm.*T1;
        [ax2,ax1,ax3] = meshgrid(1:4,1:4,1:4);
        r00 = [1 -1 0 0];%R G Y Y
        Delta = exp(-(fs-1)/fs/0.3);
        Ey1 = r00(ax2)  ... % center value
              + (r00(ax1)-r00(ax2))*Delta; % the remain of step response, if value changed
        Ey2 = r00(ax3) + (r00(ax2)-r00(ax3))*Delta;
        
        rxs1 = [rx(fs:fs:end);rx(end)];    %sampling
        PB = [1;0;0;0];
        result = zeros(length(rxs1)-1,1);
        for m = 1:length(rxs1)-1
            a = [rxs1(m),y0(m)];
            yn = rxs1(m);
            yn1 = rxs1(m+1);
            Pyyxx = exp(-((yn-Ey1).^2+(yn1-Ey2).^2)/N0);
            est = sum(T2.*Pyyxx.*PB);    % 最小差错概率估计, 对各种PB情况求和
            PB = reshape((sum(est,3)),4,1);
            PB = PB./sum(PB);
            [~,result(m)] = max(PB);
            result(m) = r00(result(m));
            
        end
        err1 = abs(result-y0);
        ErrTrial1(repTrial) = sum(err1)/T;
    end
    ErrN0(repN0) = mean(ErrTrial);
    ErrN01(repN0) = mean(ErrTrial1);
end

plot(N0_0,ErrN0)
hold on
plot(N0_0,ErrN01)
hold off



function Xf1 = PowerSpectralDensity(x,B,T,f,fs)
x = reshape(x,[],1);
f = reshape(f,1,[]);

k=((min(f)-B/2)*T:(max(f)+B/2)*T)';
% Xf=(1/T*(x.*exp(-2i*pi*t'.*k/T)))/fs;
% Xf=abs(sum(Xf)).^2;

Xf = abs(fft(x)/fs/T).^2;
Xf = Xf(mod(k,length(Xf))+1);

k1=(k<=(f+B/2)*T)&(k>=(f-B/2)*T);
Xf1=k1.*Xf;
Xf1=sum(Xf1,1);
end