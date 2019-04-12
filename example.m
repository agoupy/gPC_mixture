function x = example (delta,a,b)

c=a.^2+2.*b;
x =  (c<1/3).*(a.^2+2.*b -delta) + (c>1/3).*(-a+2.*b.^2+delta);
% x =  (c<1/3).*(a.^2+2.*b -delta) + (c>1/3).*(a.^3+2.*b+delta);

end
