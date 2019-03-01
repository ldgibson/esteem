function [ec, Vc] = VWN(p, n)
A = 0.0310907; % Eh
if n == 3
    b = 3.72744;
    c = 12.9352;
    x0 = -0.10498;
elseif n == 5
    b = 13.0720;
    c = 42.7198;
    x0 = -0.409286;
end

x = (3/(4*pi*p))^(1/6);
xi =@(q) q^2 + b*q + c;
Q = sqrt(4*c - b^2);
eta = atan(Q/(2*x+b));

ec = A*(log(x^2/xi(x)) + 2*b*eta/Q - b*x0/xi(x0)*(log((x-x0)^2/xi(x)) + ...
    2*(2*x0+b)*eta/Q));
Vc = ec - (A/3)*(c*(x-x0)-b*x*x0)/(xi(x)*(x-x0));

end