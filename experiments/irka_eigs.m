close  all;
load('CD_ss');
iter = 100;
tol = 1e-4;
r = 16;
theta = 0;
warning('off', 'all');
S = gen_shifts(r);
[Ar,Br,Cr,Sr] = irka(A,B,C,S,'std', iter, tol);
[Ai,Bi,Ci,Si] =  irka(A,B,C,S,'imag', iter, tol, theta);
G = ss(A,B,C,0);
Gr = ss(Ar,Br,Cr,0);
Gi = ss(Ai,Bi,Ci,0);

figure(1);
plot(real(eig(full(A))), imag(eig(full(A))), 'o');
hold on;
plot(real(eig(Ar)), imag(eig(Ar)), '*');

legend('Full System Poles', 'Standard Reduced System Poles', ...
    'Location', 'NorthEast');
xlabel('Re(z)');
ylabel('Im(z)');
title(sprintf('Full and Standard IRKA Eigenvalues; Reduction Order r=%d', r));
grid on;

figure(2);
plot(real(eig(full(A))), imag(eig(full(A))), 'o');
hold on;
ax = gca;
ax.ColorOrderIndex = 3;
plot(real(eig(Ai)), imag(eig(Ai)), '*');
legend('Full System Poles', 'Imaginary Reduced System Poles', ...
    'Location', 'NorthEast');
xlabel('Re(z)');
ylabel('Im(z)');
title(sprintf('Full and Imaginary Shift IRKA Eigenvalues; Reduction Order r=%d', r));
grid on;

figure(3);
bodemag(G, Gr, Gi);
legend('Full System', 'Standard Reduced System', ...
    'Imaginary Shift Reduced System', 'Location', 'NorthEast');
figure(4);
bodemag(G-Gr, G-Gi);
legend('Standard Error System', 'Imaginary Shift Error System', ...
    'Location', 'NorthEast');
figs = findobj(0, 'type', 'figure');
for k=1:length(figs)
    print(figs(k), '-depsc', sprintf('file%d.eps', k))
end