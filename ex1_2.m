Time = 800;

[Red,Green,Yellow] = meshgrid(0:Time,0:Time,0:Time);
allow = Red+Green+Yellow<=Time;
Red(~allow) = 0;
Green(~allow) = 0;
Prob = factorial(Time)./factorial(Red)./factorial(Green)./factorial(Yellow)./factorial(Time-Yellow-Red-Green)*8.^(-Time).*2.^(2*Red+Green);
Prob = sum(Prob,3);
Prob(Red+Green>Time)=0;

Num = factorial(Time)./factorial(Red)./factorial(Green)./factorial(Yellow)./factorial(Time-Yellow-Red-Green);
Num = sum(Num,3);
Num(Red(:,:,1)+Green(:,:,1)>Time) = 0;

M0 = reshape(0:80,1,1,[]);
Prob1 = repmat(Red(:,:,1)*2+Green(:,:,1),1,1,81)>repmat(M0,(Time+1),(Time+1),1);
Prob1 = Prob1.*Prob;
ProbSum = sum(sum(Prob1));
ProbSum = reshape(ProbSum,[],1);

M = find(ProbSum>0.9,1,'last');

N0 = sum(sum(Num(Red(:,:,1)*2+Green(:,:,1)>=Time)));
Nr = (ProbSum(M)-0.9)/2^(M-1)*8^Time;
N=N0+Nr;
B = log2(N)/Time


