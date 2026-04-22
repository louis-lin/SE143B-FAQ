%% wing_frequency.m
% X = span, Y = chord, Z = flap
% Section: airfoil skin + front/rear spar webs; J from inter-spar box cell only

clc; clear; close all;

%% Material
E     = 10e6;
NU    = 0.33;
G     = E / (2*(1+NU));
rho_w = 0.1;
rho   = rho_w / 386.4;

%% Wing
CHORD  = 20.0;
SPAN   = 84.0;
SKIN_T = 0.040;
SPAR_T = 0.060;
YSPAR1 = 5.0;    % front spar: 25% chord
YSPAR2 = 14.0;   % rear spar:  70% chord
MODES  = 3;

[A, Ixx, Iyy, J, Ip] = section_props(CHORD, SKIN_T, SPAR_T, YSPAR1, YSPAR2, true);

fprintf('Section (chord=%.0f", skin_t=%.3f", spar_t=%.3f")\n', CHORD, SKIN_T, SPAR_T);
fprintf('  A=%.4f in²  Ixx=%.4f in⁴  Iyy=%.4f in⁴\n', A, Ixx, Iyy);
fprintf('  J=%.4f in⁴  Ip=%.4f in⁴\n\n', J, Ip);

f_flap  = bending_modes(E, rho, SPAN, A, Ixx, MODES);
f_chord = bending_modes(E, rho, SPAN, A, Iyy, MODES);
f_tor   = torsion_modes(G, rho, SPAN, J,  Ip,  MODES);

fprintf('  Mode  |  Flap [Hz]  |  Chord [Hz]  |  Torsion [Hz]\n');
fprintf('  ------+-------------+--------------+--------------\n');
for n = 1:MODES
    fprintf('    %d   |   %7.2f   |   %8.2f   |   %8.2f\n', ...
            n, f_flap(n), f_chord(n), f_tor(n));
end

%% ═══════════════════════════════════════════════════════════════════════

function f = bending_modes(E, rho, L, A, I, N)
    bL = [1.87510, 4.69409, 7.85476, 10.99554, 14.13717];
    if N > numel(bL), bL = [bL, ((numel(bL)+1:N)-0.5)*pi]; end
    f = (bL(1:N).^2 / L^2) * sqrt(E*I / (rho*A)) / (2*pi);
end

function f = torsion_modes(G, rho, L, J, Ip, N)
    n = 1:N;
    f = ((2*n-1)*pi / (2*L)) * sqrt(G*J / (rho*Ip)) / (2*pi);
end

function [A, Ixx, Iyy, J, Ip] = section_props(chord, skin_t, spar_t, yspar1, yspar2, closed_skin)

    % --- Airfoil coords (Y=chord, Z=thickness) ---
    upper = chord * [
        0.00000,0.00000; 0.00066,0.00581; 0.00296,0.01331; 0.00688,0.02129;
        0.01236,0.02958; 0.01937,0.03807; 0.02789,0.04665; 0.03787,0.05521;
        0.04928,0.06366; 0.06206,0.07191; 0.07618,0.07987; 0.09157,0.08744;
        0.10817,0.09454; 0.12591,0.10105; 0.14478,0.10683; 0.16479,0.11180;
        0.18594,0.11591; 0.20823,0.11911; 0.23166,0.12138; 0.25623,0.12274;
        0.28190,0.12325; 0.30860,0.12294; 0.33627,0.12186; 0.36482,0.12005;
        0.39414,0.11758; 0.42415,0.11449; 0.45472,0.11085; 0.48574,0.10671;
        0.51708,0.10214; 0.54861,0.09722; 0.58019,0.09200; 0.61167,0.08655;
        0.64290,0.08094; 0.67373,0.07524; 0.70399,0.06950; 0.73353,0.06378;
        0.76219,0.05812; 0.78980,0.05258; 0.81622,0.04720; 0.84128,0.04199;
        0.86484,0.03699; 0.88676,0.03221; 0.90691,0.02765; 0.92517,0.02331;
        0.94141,0.01915; 0.95562,0.01503; 0.96797,0.01095; 0.97859,0.00715;
        0.98742,0.00396; 0.99418,0.00167; 0.99850,0.00089; 1.00000,0.00000];

    lower = chord * [
        0.00000, 0.00000; 0.00002,-0.00103; 0.00126,-0.00675; 0.00507,-0.01165;
        0.01154,-0.01641; 0.02015,-0.02102; 0.03077,-0.02538; 0.04334,-0.02937;
        0.05789,-0.03294; 0.07441,-0.03613; 0.09284,-0.03898; 0.11309,-0.04152;
        0.13505,-0.04378; 0.15860,-0.04575; 0.18361,-0.04745; 0.20997,-0.04887;
        0.23753,-0.05000; 0.26617,-0.05084; 0.29574,-0.05139; 0.32610,-0.05162;
        0.35711,-0.05153; 0.38861,-0.05110; 0.42045,-0.05031; 0.45249,-0.04910;
        0.48460,-0.04739; 0.51674,-0.04516; 0.54884,-0.04239; 0.58083,-0.03912;
        0.61264,-0.03537; 0.64419,-0.03116; 0.67547,-0.02648; 0.70000,-0.02400;
        0.72000,-0.02200; 0.74000,-0.02000; 0.76800,-0.01800; 0.79400,-0.01600;
        0.82200,-0.01400; 0.85000,-0.01200; 0.88000,-0.01000; 0.89800,-0.00800;
        0.92800,-0.00600; 0.95200,-0.00400; 0.98000,-0.00200; 1.00000, 0.00000];

    % Closed profile: upper LE→TE, lower reversed TE→LE (shared endpoints excluded)
    yp = [upper(:,1); flipud(lower(2:end-1,1))];
    zp = [upper(:,2); flipud(lower(2:end-1,2))];

    % --- Skin segments ---
    ds_sk = sqrt(diff(yp).^2 + diff(zp).^2);
    ym_sk = 0.5*(yp(1:end-1) + yp(2:end));
    zm_sk = 0.5*(zp(1:end-1) + zp(2:end));

    % --- Spar web segments (single segment each, midpoint + height) ---
    zs1t = interp1(upper(:,1), upper(:,2), yspar1);
    zs1b = interp1(lower(:,1), lower(:,2), yspar1);
    zs2t = interp1(upper(:,1), upper(:,2), yspar2);
    zs2b = interp1(lower(:,1), lower(:,2), yspar2);

    ym_sp = [yspar1;           yspar2          ];
    zm_sp = [(zs1t+zs1b)/2;   (zs2t+zs2b)/2  ];
    ds_sp = [zs1t-zs1b;        zs2t-zs2b      ];

    % --- Combined: A, centroid, Ixx, Iyy ---
    ym = [ym_sk; ym_sp];
    zm = [zm_sk; zm_sp];
    ds = [ds_sk; ds_sp];
    t  = [repmat(skin_t, numel(ds_sk), 1); spar_t; spar_t];

    dA    = t .* ds;
    A     = sum(dA);
    y_bar = sum(ym .* dA) / A;
    z_bar = sum(zm .* dA) / A;
    Ixx   = sum(t .* ds .* (zm - z_bar).^2);
    Iyy   = sum(t .* ds .* (ym - y_bar).^2);
    Ip    = Ixx + Iyy;

    % --- Bredt-Batho J: inter-spar box cell only ---
    % Box perimeter integral sum(ds/t): upper skin + lower skin + two spar webs
    in_box = ym_sk >= yspar1 & ym_sk <= yspar2;
    upper_box = in_box & zm_sk > 0;
    lower_box = in_box & zm_sk < 0;

    ds_over_t = (sum(ds_sk(upper_box)) + sum(ds_sk(lower_box))) / skin_t + ...
                (ds_sp(1) + ds_sp(2)) / spar_t;

    % Enclosed area of box cell (shoelace)
    in_u = upper(:,1) > yspar1 & upper(:,1) < yspar2;
    in_l = lower(:,1) > yspar1 & lower(:,1) < yspar2;
    ybox = [yspar1; upper(in_u,1); yspar2;          yspar2; flipud(lower(in_l,1)); yspar1];
    zbox = [zs1t;   upper(in_u,2); zs2t;   zs2b;           flipud(lower(in_l,2)); zs1b  ];
    ybox(end+1) = ybox(1); zbox(end+1) = zbox(1);  % close
    A_box = 0.5 * abs(sum(ybox(1:end-1).*zbox(2:end) - ybox(2:end).*zbox(1:end-1)));

    if closed_skin
        A_enc_full = 0.5 * abs(sum(yp(1:end-1).*zp(2:end) - yp(2:end).*zp(1:end-1)));
        J = 4 * A_enc_full^2 * skin_t / sum(ds_sk);
    else
        J = 4 * A_box^2 / ds_over_t;
    end

end
