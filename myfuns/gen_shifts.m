function [ S ] = gen_shifts( r, varargin )
%GEN_SHIFTS Summary of this function goes here
%   Detailed explanation goes here

num_varargs = length(varargin);
if num_varargs > 5
    error('inv_phi:TooManyInputs', ...
        'Specify a, b, c, d, seed1, seed2 for real in [a,b] imag in i[c,d]');
end

opt_args = {9161989, 1e-1, 1e3, 1e0, 1e5};
opt_args(1:num_varargs) = varargin;
[seed, a, b, c, d] = opt_args{:};



if mod(r, 2) == 0
    rng(seed, 'twister');
    real_S = a + (b - a).*rand(r/2, 1);
    imag_S = c + (d - c).*rand(r/2, 1);
    S = [real_S - 1i*imag_S; real_S + 1i*imag_S];
else
    
    real_S = a + (b - a).*rand((r-1)/2,1);
    real_s = a + (b - a).*rand(1, 1);
    imag_S = c + (d - c).*rand((r-1)/2, 1);
    S = [real_S - 1i*imag_S; real_S + 1i*imag_S; real_s];
end
S = sort(S, 'ascend');
end