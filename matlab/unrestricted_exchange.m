function uhf_K = unrestricted_exchange(K)
uhf_K = K;
M = length(K);
for i = 1:M
    if mod(i, 2) == 0
        % even
        uhf_K(i,1:2:M) = 0;
    else
        % odd
        uhf_K(i, 2:2:M) = 0;
    end
end

end