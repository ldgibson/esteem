function Sp = overlap_primitive(a, b, alpha, beta, A, B)

p = alpha + beta;
K_AB = exp(-alpha*beta/p*(norm(A - B))^2);
P = (alpha*A + beta*B)/p;

Sp = (pi/p)^(3/2) * K_AB * I_w(a(1), b(1), A(1), B(1), P(1), p) ...
    * I_w(a(2), b(2), A(2), B(2), P(2), p) ...
    * I_w(a(3), b(3), A(3), B(3), P(3), p);


function int1D = I_w(a_w, b_w, A_w, B_w, P_w, p)

int1D = 0;
for i = 0:(a_w+b_w)/2
    int1D = int1D + f(a_w, b_w, P_w, A_w, B_w, i) * fact2(2*i-1)/(2*p)^i;
end


function out = f(a_w, b_w, P_w, A_w, B_w, i)
out = 0;
for j = max(0, 2*i - a_w):min(2*i, b_w)
    out = out + nchoosek(a_w, 2*i - j) * nchoosek(b_w, j) ...
        * (P_w - A_w)^(a_w-2*i+j) * (P_w - B_w)^(b_w-j);
end