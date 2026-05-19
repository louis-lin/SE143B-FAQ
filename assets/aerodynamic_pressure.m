%% Wing Aerodynamic Pressure Distribution
% X = spanwise [0,l], Y = chordwise [0,c], Z = out-of-plane
% Pressure split at quarter-chord: leading [0, c/4], trailing [c/4, c]

clc; clear; clf;

%% 1. Parameters

% Full-span inputs
L_total  = 4000;   % total lift, both wings (lb)
Mx_total = 10000;  % total pitching moment, both wings (lb-in)
b        = 168;    % full wingspan (in)
c        = 20;     % chord (in)
lambda   = 1;      % taper ratio = C_tip/C_root; 1 = rectangular

% One-wing quantities
L_wing  = L_total  / 2;
Mx_star = Mx_total / 2;
l       = b / 2;

% Spanwise moment per unit span: lambda=1 collapses mx = Mx*/(4l)
mx = Mx_star / (4 * l);

%% 2. Spanwise Load Distribution Pz(x)

Pz = @(x) (L_wing/(2*l)) * ((2/pi + 1/2) - (1/pi)*(x/l).^2 - (1/(4*pi))*(x/l).^4);

figure(1); clf;
plot(0:0.1:l, Pz(0:0.1:l), 'b', 'LineWidth', 1.5);
xlabel('Spanwise X (in)'); ylabel('Pz (lb/in)');
title('Spanwise Load Distribution Pz(x)');
xlim([0 l]); grid on;

%% 3. Chordwise Pressure Fields

X_span  = 0:1:l;
Y_lead  = 0:1:c/4;
Y_trail = c/4:1:c;

[Xg_l, Yg_l] = meshgrid(X_span, Y_lead);
P_lead  =  (4/c)     * Pz(Xg_l) .* (1 - 3*(Yg_l/c)) ...
         + (24*mx/c^2)           .* (1 - 5*(Yg_l/c));

[Xg_t, Yg_t] = meshgrid(X_span, Y_trail);
P_trail =  (4/(3*c)) * Pz(Xg_t) .* (1 - Yg_t/c) ...
         - (8*mx/c^2)            .* (1 - Yg_t/c);

%% 4. Total Force (compare against RF3 at root in Abaqus)

F_lead  = trapz(X_span, trapz(Y_lead,  P_lead,  1), 2);
F_trail = trapz(X_span, trapz(Y_trail, P_trail, 1), 2);
F_total = F_lead + F_trail;

fprintf('Leading  region force: %.2f lb\n', F_lead);
fprintf('Trailing region force: %.2f lb\n', F_trail);
fprintf('Total reaction force:  %.2f lb\n', F_total);

%% 5. 3D Visualization

us = [0 0.000660 0.002960 0.006880 0.012360 0.019370 0.027890 0.037870 ...
      0.049280 0.062060 0.076180 0.091570 0.108170 0.125910 0.144780 ...
      0.164790 0.185940 0.208230 0.231660 0.256230 0.281900 0.308600 ...
      0.336270 0.364820 0.394140 0.424150 0.454720 0.485740 0.517080 ...
      0.548610 0.580190 0.611670 0.642900 0.673730 0.703990 0.733530 ...
      0.762190 0.789800 0.816220 0.841280 0.864840 0.886760 0.906910 ...
      0.925170 0.941410 0.955620 0.967970 0.978590 0.987420 0.994180 ...
      0.998500 1.000000;
      0 0.005810 0.013310 0.021290 0.029580 0.038070 0.046650 0.055210 ...
      0.063660 0.071910 0.079870 0.087440 0.094540 0.101050 0.106830 ...
      0.111800 0.115910 0.119110 0.121380 0.122740 0.123250 0.122940 ...
      0.121860 0.120050 0.117580 0.114490 0.110850 0.106710 0.102140 ...
      0.097220 0.092000 0.086550 0.080940 0.075240 0.069500 0.063780 ...
      0.058120 0.052580 0.047200 0.041990 0.036990 0.032210 0.027650 ...
      0.023310 0.019150 0.015030 0.010950 0.007150 0.003960 0.001670 ...
      0.000890 0.000100]';

ls = [0 0.000020 0.001260 0.005070 0.011540 0.020150 0.030770 0.043340 ...
      0.057890 0.074410 0.092840 0.113090 0.135050 0.158600 0.183610 ...
      0.209970 0.237530 0.266170 0.295740 0.326100 0.357110 0.388610 ...
      0.420450 0.452490 0.484600 0.516740 0.548840 0.580830 0.612640 ...
      0.644190 0.675470 0.700000 0.720000 0.740000 0.768000 0.794000 ...
      0.822000 0.850000 0.880000 0.898000 0.928000 0.952000 0.980000 ...
      1.000000;
      0 -0.001030 -0.006750 -0.011650 -0.016410 -0.021020 -0.025380 ...
     -0.029370 -0.032940 -0.036130 -0.038980 -0.041520 -0.043780 ...
     -0.045750 -0.047450 -0.048870 -0.050000 -0.050840 -0.051390 ...
     -0.051620 -0.051530 -0.051100 -0.050310 -0.049100 -0.047390 ...
     -0.045160 -0.042390 -0.039120 -0.035370 -0.031160 -0.026480 ...
     -0.024000 -0.022000 -0.020000 -0.018000 -0.016000 -0.014000 ...
     -0.012000 -0.010000 -0.008000 -0.006000 -0.004000 -0.002000 ...
     -0.000100]';

profile = [us; flip(ls,1)] * c;

figure(2); clf; hold on;
[X_surf, ~] = meshgrid(X_span, profile(:,1));
Y_surf = repmat(profile(:,1), 1, length(X_span));
Z_surf = repmat(profile(:,2), 1, length(X_span));
surf(X_surf, Y_surf, Z_surf, 'EdgeColor','none', 'FaceColor',[0.7 0.7 0.7], 'FaceAlpha',0.5);

sc = 10;  % pressure scale factor for visibility
mesh(Xg_l, Yg_l, sc*P_lead,  'EdgeColor','b', 'FaceColor','none');
mesh(Xg_t, Yg_t, sc*P_trail, 'EdgeColor','r', 'FaceColor','none');

xlabel('Spanwise X (in)'); ylabel('Chordwise Y (in)'); zlabel('Pressure (lb/in^2)');
title(sprintf('Pressure Distribution  |  Total Force = %.1f lb', F_total));
legend('Wing surface','P_{lead}','P_{trail}','Location','best');
view(-30, 30); grid on;
axis equal

%% 6. Abaqus ExpressionField Strings

abaqus_lead = sprintf(...
    '(%.4f) * (%.4f * (%.5f - %.5f*pow((X/%.0f),2) - %.5f*pow((X/%.0f),4))) * (1 - 3*(Y/%.0f)) + (%.6f) * (1 - 5*(Y/%.0f))',...
    4/c, L_wing/(2*l), (2/pi+1/2), 1/pi, l, 1/(4*pi), l, c, 24*mx/c^2, c);

abaqus_trail = sprintf(...
    '(%.5f) * (%.4f * (%.5f - %.5f*pow((X/%.0f),2) - %.5f*pow((X/%.0f),4))) * (1 - Y/%.0f) - (%.6f) * (1 - Y/%.0f)',...
    4/(3*c), L_wing/(2*l), (2/pi+1/2), 1/pi, l, 1/(4*pi), l, c, 8*mx/c^2, c);

fprintf('\nLeading edge:\n  %s\n\n', abaqus_lead);
fprintf('Trailing edge:\n  %s\n',   abaqus_trail);
