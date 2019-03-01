function [ex, Vx] = Slater(p)
% Slater exchange functional
Cx = (3/4)*(3/pi)^(1/3);
ex = -Cx*p^(1/3);
Vx = -(4/3)*Cx*p^(1/3);
end