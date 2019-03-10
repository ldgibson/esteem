function Tp = kinenergy_primitive(a, b, alpha, beta, A, B)

Tp = I_T_w(1, a, b, alpha, beta, A, B) + ... % x
    I_T_w(2, a, b, alpha, beta, A, B) + ... % y
    I_T_w(3, a, b, alpha, beta, A, B); % z

end


function I = I_T_w(w, a, b, alpha, beta, A, B)

two = [2 0 0; 0 2 0; 0 0 2];

I = beta*(2*b(w) + 1)*overlap_primitive(a, b, alpha, beta, A, B) ...
    - 2*beta^2*overlap_primitive(a, b + two(w, :), alpha, beta, A, B);

if all(b - two(w) >= 0)
    I = I - (1/2)*b(w)*(b(w)-1) * ...
        overlap_primitive(a, b - two(w, :), alpha, beta, A, B);
end    

end