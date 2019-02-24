function [J, K] = eerepulsion(ERI, P)
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
                K_ = K_ + 0.5*P(kap,lam)*ERI(mu, kap, lam, nu);
            end
        end
        J(mu,nu) = J_;
        K(mu,nu) = K_;
        if mu ~= nu
            J(nu,mu) = J_;
            K(nu,mu) = K_;
        end
    end
end
% Vee = J - K;
end