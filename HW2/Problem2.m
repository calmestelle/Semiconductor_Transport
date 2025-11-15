%% Problem 2
% *Using the tight-binding (LCAO) approach, compute the energy band formed by 
% atomic s-orbitals in an fcc lattice with a monoatomic basis. Plot the energy 
% vs q along L -* $\Gamma$ *- X - U -* $\Gamma$

clc
clear
close all

gamma1 = -1;

Gamma = [0, 0, 0];
X = [2*pi, 0, 0];
L = [pi, pi, pi];
U = [2*pi, pi/2, pi/2];
% Since a cancels out when calculating the energy, plus we don't have a
% value to set it to, I'm going to leave it out of this solution.

k_LG = [linspace(L(1), Gamma(1), 100); linspace(L(2), Gamma(2), 100); linspace(L(3), Gamma(3), 100)];
k_GX = [linspace(Gamma(1), X(1), 100); linspace(Gamma(2), X(2), 100); linspace(Gamma(3), X(3), 100)];
k_XU = [linspace(X(1), U(1), 100); linspace(X(2), U(2), 100); linspace(X(3), U(3), 100)];
k_UG = [linspace(U(1), Gamma(1), 100); linspace(U(2), Gamma(2), 100); linspace(U(3), Gamma(3), 100)];

E_LG = 2 * gamma1 * (cos((k_LG(1,:) + k_LG(2,:))/2) + cos((k_LG(1,:) - k_LG(2,:))/2) + ...
    (cos((k_LG(1,:) + k_LG(3,:))/2) + cos((k_LG(1,:) - k_LG(3,:))/2)) + ...
    (cos((k_LG(2,:) + k_LG(3,:))/2) + cos((k_LG(2,:) - k_LG(3,:))/2)));
E_GX = 2 * gamma1 * (cos((k_GX(1,:) + k_GX(2,:))/2) + cos((k_GX(1,:) - k_GX(2,:))/2) + ...
    (cos((k_GX(1,:) + k_GX(3,:))/2) + cos((k_GX(1,:) - k_GX(3,:))/2)) + ...
    (cos((k_GX(2,:) + k_GX(3,:))/2) + cos((k_GX(2,:) - k_GX(3,:))/2)));
E_XU = 2 * gamma1 * (cos((k_XU(1,:) + k_XU(2,:))/2) + cos((k_XU(1,:) - k_XU(2,:))/2) + ...
    (cos((k_XU(1,:) + k_XU(3,:))/2) + cos((k_XU(1,:) - k_XU(3,:))/2)) + ...
    (cos((k_XU(2,:) + k_XU(3,:))/2) + cos((k_XU(2,:) - k_XU(3,:))/2)));
E_UG = 2 * gamma1 * (cos((k_UG(1,:) + k_UG(2,:))/2) + cos((k_UG(1,:) - k_UG(2,:))/2) + ...
    (cos((k_UG(1,:) + k_UG(3,:))/2) + cos((k_UG(1,:) - k_UG(3,:))/2)) + ...
    (cos((k_UG(2,:) + k_UG(3,:))/2) + cos((k_UG(2,:) - k_UG(3,:))/2)));

E = [E_LG; E_GX; E_XU; E_UG];
E = E(:);

figure;
hold on;
grid on;
% L -> 
subplot(2, 2, 1);
plot(E_LG, 'b', 'LineWidth', 2);
title('Energy Band vs q: L -> ');
ylabel('Energy E(q)');

%  -> X
subplot(2, 2, 2);
plot(E_GX, 'g', 'LineWidth', 2);
title('Energy Band vs q:   X');
ylabel('Energy E(q)');

% X -> U
subplot(2, 2, 3);
plot(E_XU, 'r', 'LineWidth', 2);
title('Energy Band vs q: X -> U');
xlabel("q (Reciprocal space)");
ylabel('Energy E(q)');

% U -> 
subplot(2, 2, 4);
plot(E_UG, 'r', 'LineWidth', 2);
title('Energy Band vs q: U -> ');
xlabel("q (Reciprocal space)");
ylabel('Energy E(q)');
hold off