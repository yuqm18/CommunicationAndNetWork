
a = 1;%重建电平 -a,0,a
b = 0.5;%量化电平 -b，b

err = @(a,b) (1+a^2*2*(1-normcdf(b))-4*a/sqrt(2*pi)*exp(-b^2/2));

% b = 0.5, find best a
besta = @(b)2/sqrt(2*pi).*exp(-b.^2/2)/( (normcdf(b)-0.5)*(-2)+1);
E1 = err(besta(0.5),0.5);

% find best a,b
err1 = @(x) err(x(1),x(2));
x = fminsearch(err1,[0.5,1]);
err1(x)