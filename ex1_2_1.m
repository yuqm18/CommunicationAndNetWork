Time = 40;

[Red,Green] = meshgrid(0:Time,0:Time);
allow = Red+Green<=Time;
Red(~allow) = 0;
Green(~allow) = 0;
Prob = factorial(Time)./factorial(Red)./factorial(Green)./factorial(Time-Red-Green).*2.^(Time-Red-Green).*8.^(-Time).*2.^(2*Red+Green);
Prob(~allow)=0;

Num = factorial(Time)./factorial(Red)./factorial(Green)./factorial(Time-Red-Green).*2.^(Time-Red-Green);
Num(~allow) = 0;

M0 = reshape(0:Time*2,1,1,[]);
Prob1 = repmat(Red(:,:,1)*2+Green(:,:,1),1,1,Time*2+1)>repmat(M0,(Time+1),(Time+1),1);
Prob1 = Prob1.*Prob;
ProbSum = sum(sum(Prob1));
ProbSum = reshape(ProbSum,[],1);

M = find(ProbSum>0.9,1,'last');

N0 = sum(sum(Num(Red*2+Green>=M)));
Nr = (ProbSum(M)-0.9)/2^(M-1)*8^Time;
N=N0+Nr;
B = log2(N)/Time


