function [ iS ] = inv_phi( S, varargin )
%INV_PHI Summary of this function goes here
%   Detailed explanation goes here
num_varargs = length(varargin);
if num_varargs > 1
    error('inv_phi:TooManyInputs', ...
        'Only specify inverse angle in radians');
end

opt_args = {-pi};
opt_args(1:num_varargs) = varargin;
theta = opt_args{:};
S = full(S);
numer = poly(S);
denom = poly(-S);

s = exp(1i * theta);
iS = roots(numer - s * denom);
iS = 0 + 1i*imag(iS);

iS(abs(imag(iS)) < 1e-8) = 0;
iS_pos = sort(abs(iS(imag(iS) > 0)), 'ascend');
iS_neg = sort(abs(iS(imag(iS) < 0)), 'ascend');
sprintf('pos length:\t%d\nneg length:\t%d', length(iS_pos), length(iS_neg));
iS = sort(iS, 'ascend');
save('iS.mat', 'iS');
iS_zero = iS(imag(iS) == 0);
iS_avg = (iS_pos + iS_neg) / 2;
iS = 1i * [iS_avg; -iS_avg; iS_zero];
iS = sort(iS, 'ascend');

if length(iS) < length(S)
    iS = [iS; 0];
    iS = sort(iS, 'ascend');

end