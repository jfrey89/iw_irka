function [ v, w ] = krylov_zero( A, B, C, n_zeros )
%KRYLOV_ZERO Summary of this function goes here
%   Detailed explanation goes here
v = -A\B;
w = -A.'\C.';
if n_zeros == 0
    return
else
    for i = 1:n_zeros
        v = -A\v;
        w = -A.'\w;
    end
end

end