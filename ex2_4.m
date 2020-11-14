syms  p e
q = 1-p;
f = 1-e;

pd = [q^2,p*q,p*q,p*p].';

py = [f*f,e*f,e*f,e*e;
      e*e,e*f,e*f,f*f;
      e*f,e*e,f*f,e*f;
      e*f,f*f,e*e,e*f];%第i行，给定d=i，p(Y|D)
  
py1 = sum(pd.*py);%p(Y)

Hy = -sum(py1.*log2(py1));

pyd = (pd.*py);%P(YD)
pyd1(1,:) = sum(pyd(1:2,:),1);
pyd1(2,:) = sum(pyd(3:4,:),1);%P(YD1)


Hyd1 = -sum([q;p].*sum(pyd1.*log2(pyd1),2));

ID1Y  = Hy-Hyd1;

IDIY1 = sum(sum(pyd1.*(log2(pyd1)-log2([q;p])-log2(py1))));

simplify(ID1Y-IDIY1)
simplify(ID1Y+IDIY1)

    

  
  