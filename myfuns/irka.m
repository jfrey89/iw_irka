function [ Ar, Br, Cr, S_iter] = irka( A, B, C, S, varargin)
%IRKA Summary of this function goes here
%   Detailed explanation goes here
% warning('off', 'all');
num_varargs = length(varargin);
if num_varargs > 4
    error('irka:TooManyInputs', ...
        'Either std or imag');
end

% set default for optional inputs
opt_args = {'std', 100, 1e-4, -pi};

% now put these defaults into the values_to_use cell array,
% and overwrite the ones specified in varargin
opt_args(1:num_varargs) = varargin;

% place optional args in memorable variable names
[shift_type, max_iter, tol, theta] = opt_args{:};

iter = 1;
r = length(S);
S_iter = zeros(r, max_iter + 1);
shift_error = 1;

if strcmp(shift_type, 'std')
%     fprintf('\nSTANDARD SHIFTS\n');
    S_iter(:, iter) = sort(S, 'ascend');
    [V, W] = krylov_red(A, B, C, S);
    [Ar, Br, Cr, S] = model_red(A, B, C, V, W);
    x_old = vectorizeS(S);
    
    while shift_error > tol && iter < max_iter
        iter = iter + 1;
        S_iter(:, iter) = sort(S, 'ascend');
        [V, W] = krylov_red(A, B, C, S);
        [Ar, Br, Cr, S] = model_red(A, B, C, V, W);
        x_new = vectorizeS(S);
        shift_error = norm(x_new - x_old) / norm(x_new);
        if (mod(iter, 15) == 0)
%             fprintf('Shift Error:\t%f\n', shift_error);
        end
        x_old = x_new;
    end
    
    S_iter = S_iter(:, 1:iter);
    
elseif strcmp(shift_type, 'imag')
%     fprintf('\nIMAGINARY SHIFTS\n');
    S_iter(:, iter) = S;
    S = inv_phi(S, theta);
    [V, W] = krylov_red(A, B, C, S);
    [Ar, Br, Cr, S] = model_red(A, B, C, V, W);
    x_old = vectorizeS(S);
    
    while shift_error > tol && iter < max_iter
        iter = iter + 1;
        S_iter(:, iter) = S;
        S = inv_phi(S, theta);
        [V, W] = krylov_red(A, B, C, S);
        [Ar, Br, Cr, S] = model_red(A, B, C, V, W);
        x_new = vectorizeS(S);
        shift_error = norm(x_new - x_old) / norm(x_new);
        if (mod(iter, 15) == 0)
%             fprintf('Shift Error:\t%f\n', shift_error);
        end
        x_old = x_new;
    end
    S_iter = S_iter(:, 1:iter);
else
    error('irka:ShiftsNotSpecified', ...
        'Specify std for classical IRKA or imag for imshiftIRKA');
end
% warning('on', 'all');
end