function [ vS ] = vectorizeS( S )
%VECTORIZES Summary of this function goes here
%   Detailed explanation goes here
r = length(S);
S = sort(S, 'ascend');
vS = zeros(r, 1);

i = 1;
while i <= r
    if abs(imag(S(i)) / abs(S(i))) < 1e-5
        vS(i) = abs(S(i));
        i = i + 1;
    else
        vS(i) = real(S(i));
        vS(i + 1) = abs(imag(S(i)));
        i = i + 2;
    end
end

end

