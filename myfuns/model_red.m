function [ Ar, Br, Cr, S ] = model_red( A, B, C, V, W )
%MODEL_REDUCE Summary of this function goes here
%   Detailed explanation goes here
Ar = W.' * A * V;
Br = W.' * B;
Cr = C * V;
S = -eig(full(Ar));
S = sort(S, 'ascend');
end