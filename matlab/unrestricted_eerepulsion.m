function [J, K] = unrestricted_eerepulsion(ERI, P)
M = length(ERI);
J = zeros(M);
K = zeros(M);

for mu = 1:M
    for nu = mu:M
        J_ = 0;
        K_ = 0;
        for kap = 1:M
            for lam = 1:M
                J_ = J_ + P(kap,lam)*ERI(mu, nu, lam, kap);
                K_ = K_ + P(kap,lam)*ERI(mu, kap, lam, nu);
            end
        end
        J(mu,nu) = J_;
        if mod(mu, 2) == 0 && mod(nu, 2) == 0
            K(mu,nu) = K_;
        elseif mod(mu-1, 2) == 0 && mod(nu-1, 2) == 0
            K(mu,nu) = K_;
        else
            K(mu,nu) = 0;
            K_ = 0;
        end
        if mu ~= nu
            J(nu,mu) = J_;
            K(nu,mu) = K_;
        end
    end
end
% Vee = J - K;
end