close all;
load('is');
tol = 1e-4;

r1 = 8;
r2 = 10;
sysf = ss(A,B,C,0);
a = 1e-1; b = 1e3; c = 1e0; d = 1e5;
N = 15;
errs = zeros(N, 2);
norm_sysf = norm(sysf);

for iter=1:N
    passed1 = 0;
    passed2 = 0;
    
    while passed1 < 1
        S1 = gen_shifts(r1, a, b, c, d);
        
        [Ar, Br, Cr, Sr] = irka(A, B, C, S1, 'imag', iter, tol);
        sysr = ss(Ar, Br, Cr, 0);
        err = norm(sysf - sysr)/norm_sysf;
        
        if isinf(err) < 1
            passed1 = 1;
            errs(iter, 1) = err;
        end
    end
    while passed2 < 1
        S2 = gen_shifts(r2, a, b, c, d);
        [Ar, Br, Cr, Sr] = irka(A, B, C, S2, 'imag', iter, tol);
        sysr = ss(Ar, Br, Cr, 0);
        err = norm(sysf - sysr)/norm_sysf;
        
        if isinf(err) < 1
            passed2 = 1;
            errs(iter, 2) = err;
        end
    end
end
figure1 = figure(1); clf;
plot(errs(:,1), '-^');
hold on;
plot(errs(:,2), '-s');
ylabel_str = '$\frac{||G-G_r||_2}{||G||_2}$';
ylabel(ylabel_str, 'Interpreter', 'LaTeX');
xlabel('# iterations');
legend('r = 8', 'r = 10', 'Location', 'NorthEast');
print -depsc2 error_vs_iter.eps