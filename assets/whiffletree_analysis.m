%% Whiffletree Two-Point Load Analysis
%  Approximates Schrenk spanwise distribution with two point loads.
%  Matching conditions: total force + root bending moment (exact at root).
%  Students edit STUDENT INPUTS block only.

%% Fixed Parameters
L_total = 4000;          % total lift, both wings (lb)
b       = 168;           % full wingspan (in)
l       = b / 2;         % half-span (in)
L_wing  = L_total / 2;  % one-wing lift (lb)

Pz = @(x) (L_wing/(2*l)) * ((2/pi + 0.5) - (x/l).^2/pi - (x/l).^4/(4*pi));

%% STUDENT INPUTS
x1 = 28;   % inboard load location (in from root)
x2 = 63;   % outboard load location (in from root)
% Spacing = x2 - x1

%% Solve for P1, P2
xv  = linspace(0, l, 2000);
pzv = Pz(xv);

F0 = trapz(xv, pzv);           % total force = L_wing (verify)
M0 = trapz(xv, xv .* pzv);    % first moment about root

P  = [1, 1; x1, x2] \ [F0; M0];
P1 = P(1);
P2 = P(2);

fprintf('\n--- Whiffletree Solution ---\n');
fprintf('x1 = %5.1f in  |  P1 = %7.1f lb\n', x1, P1);
fprintf('x2 = %5.1f in  |  P2 = %7.1f lb\n', x2, P2);
fprintf('P1 + P2 = %.1f lb  (target: %.1f lb)\n', P1+P2, L_wing);
fprintf('Spacing = %.1f in\n\n', x2-x1);

%% Bending Moment Diagrams
xp  = linspace(0, l, 800);
pzp = Pz(xp);

Q   = cumtrapz(xp, pzp);
Q   = Q(end) - Q;
R   = cumtrapz(xp, xp .* pzp);
R   = R(end) - R;
Md  = R - xp .* Q;

Mp  = (xp <= x1) .* P1 .* (x1 - xp) + (xp <= x2) .* P2 .* (x2 - xp);

rms_err = sqrt(mean((Mp - Md).^2)) / Md(1) * 100;
%% Plots
figure('Position', [100 80 860 680], 'Color', 'w');

% --- Load Distribution ---
subplot(2,1,1);
fill([xv, fliplr(xv)], [pzv, zeros(size(pzv))], ...
    [0.75 0.88 1], 'EdgeColor', [0 0.4 0.8], 'LineWidth', 1.5);
hold on;
ymax = max(pzv) * 1.15;
plot([x1 x1], [0 ymax], 'r--', 'LineWidth', 1.5);
plot([x2 x2], [0 ymax], 'm--', 'LineWidth', 1.5);
text(x1, ymax * 0.88, sprintf('x_1 = %g"\nP_1 = %.0f lb', x1, P1), ...
    'Color', 'r', 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
text(x2, ymax * 0.88, sprintf('x_2 = %g"\nP_2 = %.0f lb', x2, P2), ...
    'Color', [0.6 0 0.8], 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
xlabel('Spanwise Position (in)');
ylabel('p_z  (lb/in)');
title('Schrenk Load Distribution with Point Load Locations');
xlim([0 l]); ylim([0 ymax]); grid on; box on;

% --- Moment Diagram ---
subplot(2,1,2);
plot(xp, Md, 'b-',  'LineWidth', 2.5); hold on;
plot(xp, Mp, 'r--', 'LineWidth', 2);
xline(x1, ':', 'Color', 'r',       'LineWidth', 1);
xline(x2, ':', 'Color', [0.6 0 0.8], 'LineWidth', 1);
xlabel('Spanwise Position (in)');
ylabel('Bending Moment (lb-in)');
title(sprintf('Bending Moment Diagram  |  x_1=%g", x_2=%g"  |  RMS error = %.1f%%', ...
    x1, x2, rms_err));
legend('Schrenk (distributed)', 'Two-point approximation', 'Location', 'northeast');
xlim([0 l]); grid on; box on;

fprintf('Moment at root (Schrenk):    %.1f lb-in\n', Md(1));
fprintf('Moment at root (two-point):  %.1f lb-in\n', Mp(1));
fprintf('RMS moment error:            %.2f%%\n', rms_err);
