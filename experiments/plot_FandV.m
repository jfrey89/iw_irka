load CD_ss;
gr = (1 + sqrt(5))/2;
sysf = ss(A,B,C,0);
iter = 100;
tol = 1e-4;
r = 16;
a = 1e-1; b = 1e3;
c = 1e0; d = 1e5;
res = 5e4;
% S = gen_shifts(r, a, b, c, d);
% [Ar, Br, Cr] = irka(A, B, C, S, 'std', iter, tol);
sysr = ss(Ar, Br, Cr, 0);
z = eig(Ar);
y = -real(z) + 1i*imag(z);
% y = inv_phi(y);

figure(1), clf;
hold on;
axis equal
% xbnd = gr*max(abs(real(y)));
% xbnd = gr * mean(abs(real(y)));
xbnd = 500;
ybnd = xbnd;
axis([-xbnd xbnd -ybnd ybnd])
plot([0 0],[-xbnd xbnd],'m-')

figure(2), clf;
T = exp(linspace(0, 2i*pi, res));
plot(T, 'm-');
hold on;
axis equal, axis([-gr gr -gr gr]);
num = poly(y);
den = poly(z);

k = 4;
dvec = zeros(2*k+1,1);
j = 1;
for i = -k:k
    dvec(j) = gr^i;
    j = j + 1;
end
dvec = [dvec(1:k); dvec(k+2:end)];
TT = exp(linspace(0, 2i*pi, res));

for j = 1:length(dvec)
    delta = dvec(j);
    Td = delta*TT;
    Zhat = zeros(r * res, 1);
    idx = 0;
    
    for k = 1:length(Td)
        zhat = roots(num - Td(k)*den);
        Zhat(idx + 1:k*r) = zhat;
        idx = idx + r;
    end
    [~, idx] = sort(imag(Zhat));
    Zhat = Zhat(idx);
    figure(2);
    hold on;
    plot(delta*T, 'k-', 'linewidth', 3)
    
    figure(1);
    hold on;
    plot(real(Zhat), imag(Zhat), 'k.', 'markersize', 3);
end

thvec = -pi + 2*pi/r * linspace(0, r-1, r);
for j = 1:length(thvec)
    Rth = linspace(0, 100, res) * exp(1i*thvec(j));
    
    Zhat = zeros(res, r);
    for k = 1:length(Rth)
        zhat = roots(num - Rth(k)*den);
        [~, idx] = sort(imag(zhat));
        Zhat(k, :) = zhat(idx).';
    end
    figure(2);
    plot(real(Rth), imag(Rth), 'g-', ...
        'markersize', 3, 'linewidth', 3);
    figure(1);
    hold on;
    for k = 1:r
        zhat = Zhat(:, k);
        plot(real(zhat), imag(zhat), 'g.', ...
            'markersize', 3, 'linewidth', 3);
    end
end

Zhat = zeros(r*r, 1);
idx = 0;
for j = 1:length(thvec)
    pt = 1*exp(1i*thvec(j));
    zhat = roots(num - pt*den);
    if length(zhat) < r
        zhat = [zhat; min(zhat)];
        zhat = sort(zhat, 'ascend');
    end
    Zhat(idx+1:j*r) = zhat;
    idx = idx + r;
end
figure(1);
hold on;
% plot(real(Zhat), imag(Zhat), 'ko', 'markersize', 6);
% phi_y = inv_phi(y);
% plot(real(phi_y), imag(phi_y), 'o', 'markersize', 6, ...
%     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');

figure(1);
plot(z, 'o', 'markersize', 6, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
hold on;

plot(y, 'o', 'markersize', 6, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');

% for j=2:length(thvec)
%     th = thvec(j);
%     phi_y = inv_phi(y, th);
%     plot(real(phi_y), imag(phi_y), 'o', 'markersize', 5, ...
%         'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y');
% end

% figs = findobj(0, 'type', 'figure');
% for k=1:length(figs)
%     print(figs(k), '-dpng', sprintf('file%d.png', k))
% end