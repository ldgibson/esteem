function prim_int = VRR(m, a, b, alpha, beta, A, B, C)
% prim_int = VRR(m, a, b, alpha, beta, A, B, C)
%
% Input:
%   m           order of auxiliary function
%   a, b        vector of Cartesian exponents
%   alpha, beta radial exponents of primitives
%   A, B, C     vectors of Cartesian coordinates of nuclei
% Output:
%   prim_int    integral of primitives; [a|r_c^-1|b]

p = alpha + beta;
P = (alpha*A + beta*B)/p;
K_AB = exp(-alpha*beta/p*norm(A - B)^2);
T = p * norm(P - C)^2;

one = [1 0 0; 0 1 0; 0 0 1];
two = [2 0 0; 0 2 0; 0 0 2];

if all(a == 0) && all(b == 0)
    prim_int = (2*pi/p)*K_AB*boysF(m, T);
elseif any(a > 0)
    if a(1) ~= 0
        w = 1;
    elseif a(2) ~= 0
        w = 2;
    else
        w = 3;
    end

    % if a(w) >= 2 and b(w) == 1
    if all(a - two(w, :) >= 0) && all(b - one(w, :) >= 0)
        prim_int =  (P(w) - A(w))*VRR(m, a - one(w, :), b, alpha, beta, A, B, C) ...
            + (C(w) - P(w))*VRR(m+1, a - one(w, :), b, alpha, beta, A, B, C) ...
            + (a(w) - 1)/(2*p)*(VRR(m, a - two(w, :), b, alpha, beta, A, B, C) ...
            - VRR(m+1, a - two(w, :), b, alpha, beta, A, B, C)) ...
            + (b(w)/(2*p))*(VRR(m, a - one(w, :), b - one(w, :), alpha, beta, A, B, C) ...
            - VRR(m+1, a - one(w, :), b - one(w, :), alpha, beta, A, B, C));
    % if a(w) == 2
    elseif all(a - two(w, :) >= 0)
        prim_int =  (P(w) - A(w))*VRR(m, a - one(w, :), b, alpha, beta, A, B, C) ...
            + (C(w) - P(w))*VRR(m+1, a - one(w, :), b, alpha, beta, A, B, C) ...
            + (a(w) - 1)/(2*p)*(VRR(m, a - two(w, :), b, alpha, beta, A, B, C) ...
            - VRR(m+1, a - two(w, :), b, alpha, beta, A, B, C));
    % if a(w) == 1 and b(w) >= 1
    elseif all(b - one(w, :) >= 0)
        prim_int = (P(w) - A(w))*VRR(m, a - one(w, :), b, alpha, beta, A, B, C) ...
            + (C(w) - P(w))*VRR(m+1, a - one(w, :), b, alpha, beta, A, B, C) ...
            + (b(w)/(2*p))*(VRR(m, a - one(w, :), b - one(w, :), alpha, beta, A, B, C) ...
            - VRR(m+1, a - one(w, :), b - one(w, :), alpha, beta, A, B, C));
    % if a(w) == 1
    else
        prim_int = (P(w) - A(w))*VRR(m, a - one(w, :), b, alpha, beta, A, B, C) ...
            + (C(w) - P(w))*VRR(m+1, a - one(w, :), b, alpha, beta, A, B, C);
    end
elseif any(b > 0)
    if b(1) ~= 0
        w = 1;
    elseif b(2) ~= 0
        w = 2;
    else
        w = 3;
    end
    % if b(w) >= 2
    if all(b - two(w, :) >= 0)
        prim_int = (P(w) - B(w))*VRR(m, a, b - one(w, :), alpha, beta, A, B, C) ...
            + (C(w) - P(w))*VRR(m+1, a, b - one(w, :), alpha, beta, A, B, C) ...
            + (b(w) - 1)/(2*p)*(VRR(m, a, b - two(w, :), alpha, beta, A, B, C) ...
            - VRR(m+1, a, b - two(w, :), alpha, beta, A, B, C));
    % if b(w) == 1
    else
        prim_int = (P(w) - B(w))*VRR(m, a, b - one(w, :), alpha, beta, A, B, C) ...
            + (C(w) - P(w))*VRR(m+1, a, b - one(w, :), alpha, beta, A, B, C);
    end
end
end