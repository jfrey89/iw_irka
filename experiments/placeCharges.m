load iss;
G = ss(A,B,C,0);
r = 8;
S = gen_shifts(r);
[Ar,Br,Cr] = irka(A,B,C,S);
Gr = ss(Ar,Br,Cr,0);
z = eig(Ar);
y = -real(z)+1i*imag(z);
num = poly(y);
den = poly(z);


thvec = -pi + 2*pi/r * linspace(0, r-1, r);


for j = 1:length(thvec)
    Rth = 1;
    zhat = roots(num - Rth*den);
end