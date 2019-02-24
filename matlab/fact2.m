function f2 = fact2(n)
if mod(n, 2) == 0
    f2 = prod([2:2:n]);
else
    f2 = prod([1:2:n]);
end
end