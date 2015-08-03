function [ V, W ] = krylov_red( A, B, C, S )
%KRYLOV_RED Summary of this function goes here
%   Detailed explanation goes here

n = size(A, 2);
E = eye(n);
r = length(S);
V = zeros(n, r);
W = zeros(n, r);
S = sort(S, 'ascend');

i = 1;
n_zeros = 0;

while i <= r
    if abs(S(i)) == 0
        [v, w] = krylov_zero(A, B, C, n_zeros);
        V(:, i) = v;
        W(:, i) = w;
        i = i + 1;
        n_zeros = n_zeros + 1;
%     elseif abs(imag(S(i))) / abs(S(i)) <= 1e-5
%         v = (real(S(i)) * E - A)  \ B;
%         V(:, i) = real(v);
%         w = (real(S(i)) * E - A).' \ C.';
%         W(:, i) = real(w);
%         i = i + 1;
    elseif abs(imag(S(i))) / abs(S(i)) > 1e-5 && i < r
        v = (S(i) * E - A)  \ B;
        V(:, i) = real(v);
        V(:, i + 1) = imag(v);
        w = (S(i) * E - A).' \ C.';
        W(:, i) = real(w);
        W(:, i + 1) = imag(w);
        i = i + 2;
    else
        v = (real(S(i)) * E - A)  \ B;
        V(:, i) = real(v);
        w = (real(S(i)) * E - A).' \ C.';
        W(:, i) = real(w);
        i = i + 1;
    end
end

[V, ~, ~] = svd(V, 0);
[W, ~, ~] = svd(W, 0);

W = W / (W.' * V)';

end