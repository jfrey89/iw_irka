close all;
load('CD_ss');
tol = 1e-4;
r_min = 2;
r_max = 40;
iter = 100;
theta = -pi;
warning('off', 'all');

errors = zeros(length(r_min:2:r_max), 3);
tic
sysf = prescale(ss(A, B, C, 0));
norm_sysf = norm(sysf);

a = 1e-1; b = 1e3;
c = 1e0; d = 1e5;

i = 1;
for r = r_min:r_max
    if mod(r, 2) == 0
        fprintf('r = %d\n', r)
        passed1 = 0;
        passed2 = 0;
        while passed1 + passed2 < 2
            S = gen_shifts(r, randi(2^32), a, b, c, d);
            while passed1 < 1
%                 S = gen_shifts(r);
                [Ar, Br, Cr, S_iter] = irka(A, B, C, S, ...
                    'std', iter, tol);
                sysr = ss(Ar, Br, Cr, 0);
                err_r = norm(sysf - sysr);
                
                if isinf(err_r) < 1
                    fprintf('Std. Passed\n');
                    passed1 = 1;
                else
                    fprintf('Std. Failed\n');
                end
            end
            passed2 = 0;
            while passed2 < 1
%                 S = gen_shifts(r);
                [Ari, Bri, Cri, S_iteri] = irka(A, B, C, S, ...
                    'imag', iter, tol, theta);
                sysi = ss(Ari, Bri, Cri, 0);
                err_i = norm(sysf - sysi);
                if isinf(err_i) < 1
                    fprintf('Equiv. Passed\n');
                    passed2 = 1;
                else
                    fprintf('Equiv. Failed\n');
                end
            end
%             
            if passed1 + passed2 == 2
                errors(i, :) = [r err_r/norm_sysf err_i/norm_sysf];
                i = i + 1;
                
            end
        end
    end
end
toc
figure;
semilogy(errors(:, 1), errors(:, 2), '-^', ...
    errors(:, 1), errors(:, 3), '-s');
ylabel('$\frac{||G - G_r||_2}{||G||_2}$', 'Interpreter', 'LaTeX');
xlabel('r');
ax = gca;
ax.XTick = 0:4:r_max;
legend('Standard IRKA', 'Imaginary Shift IRKA', 'Location', 'NorthEast');
% set(gcf, 'PaperPositionMode', 'auto');
% figure;
% plot(real(eig(full(A))), imag(eig(full(A))), '.', ...
%     real(eig(Ar)), imag(eig(Ar)), '.', ...
%     'MarkerSize', 10);
% legend('full eigenvalues', 'standard reduced eigenvalues', ...
%     'Location', 'West');
% hold on;
% figure;
% plot(real(eig(full(A))), imag(eig(full(A))), '.', ...
%     'MarkerSize', 10);
% hold on
% ax = gca;
% ax.ColorOrderIndex = 3;
% plot(real(eig(Ari)), imag(eig(Ari)), '.', ...
%     'MarkerSize', 10);
% legend('full eigenvalues', 'imaginary shifts reduced eigenvalues', ...
%     'Location', 'West');
% figs = findobj(0, 'type', 'figure');
% for k=1:length(figs)
%     print(figs(k), '-deps', sprintf('file%d.eps', k))
% end
print '-depsc2' 'error_vs_r.eps'
warning('on', 'all');
