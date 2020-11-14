L = 100000;
data = rand(L,1)>0.25;
vol = data.*(3+1i)+(1-data).*(-1-1i);
err = ones(10,1);
for k = 1:10
    n2 = k;
vol1 = vol+random('Normal',0,n2/2,[L,1])+1i*random('Normal',0,n2/2,[L,1]);
gamma = 8/n2;

data_out = real(vol1)+imag(vol1)/2>1-1/gamma*log(3);

err(k) = min(sum(data_out~=data)/L,err(k));

err0(k) = (1-normcdf(sqrt(gamma/20)*(5-2/gamma*log(3))))*0.25+(1-normcdf(sqrt(gamma/20)*(5+2/gamma*log(3))))*0.75;
end
plot(err)
hold on
plot(err0)
hold off