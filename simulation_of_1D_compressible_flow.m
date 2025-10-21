% duct_flow.m
% =======================================================================
%  duct_flow.m
%  Steady 1-D compressible flow solver for ducts with:
%    - area change
%    - friction
%    - heat addition
%
%  Author: Mustafa Taha
%  Date: April 2024
%
%  Description:
%  ------------
%  This MATLAB script solves the governing equations of 1-D compressible flow
%  using finite-difference marching along a specified duct geometry. 
%  It models:
%     1. Isentropic (area change only)
%     2. Fanno (with friction)
%     3. Rayleigh + Fanno (with friction and heat)
%
%  Outputs: Mach number, velocity, temperature, pressure, density distributions
%           along the duct and visual plots.
%
% =======================================================================



% Steady 1-D compressible flow in a duct with area change, friction, heat
% Usage: run the script. Adjust inputs at the top.

clearvars; close all; clc;

% -------------------- INPUT --------------------
NIU = 1;          % 1 = axisymmetric, 0 = 2D (unused here but kept)
XL = 1;           % duct length (m)
imax = 10001;     % number of grid points
MIN = 0.137;      % inlet Mach (subsonic case). Use e.g. 2.2 for supersonic
TTIN = 840;       % total temperature at inlet (K)
PTIN = 3721.7e3;  % total pressure at inlet (Pa)
gama = 1.4; RGAS = 287;
CP = gama * RGAS / (gama - 1);
% ------------------ END INPUT -------------------

% Derived constants
dx = XL / (imax - 1);

% Precompute constants used later
R1 = .5 * (gama + 1) / gama;
R2 = .5 * (gama - 1) / gama;
R3 = 2 * gama * (gama - 1) / (gama + 1);

% ---------------- DUCT AREA DISTRIBUTION ----------------
% Example convergent-divergent defined piecewise by cubic interpolation
A1 = 0.5612; A2 = 0.1403; A3 = 0.5612;
x1 = 0; x2 = 0.2113; x3 = 1;

% polynomial coefficients for the two segments (Hermite-like)
c0 = A1;
c2 = 3*(A2-A1)/(x2-x1)^2;
c3 = -2*(A2-A1)/(x2-x1)^3;
d0 = A2;
d2 = 3*(A3-A2)/(x3-x2)^2;
d3 = -2*(A3-A2)/(x3-x2)^3;

% build arrays
a_x = linspace(0,XL,imax);
A = nan(1,imax);
for i = 1:imax
    X = a_x(i);
    if X <= x2
        A(i) = c0 + c2*(X-x1)^2 + c3*(X-x1)^3;
    else
        A(i) = d0 + d2*(X-x2)^2 + d3*(X-x2)^3;
    end
end
a_D = sqrt(4*A/pi);    % equivalent diameter

% ---------------- Run 3 cases ----------------
% cases:
% 1: area change only (isentropic)
% 2: area change + friction
% 3: area change + friction + heat
results = cell(3,1);

for icase = 1:3
    switch icase
        case 1
            Q_total = 0;
            FR = 0;
        case 2
            Q_total = 0;
            FR = 0.02;
        case 3
            Q_total = 4e5;
            FR = 0.02;
    end

    dQ = Q_total * dx / XL;  % heat added per cell (J/kg per step approximation)

    % initial conditions (filled to arrays)
    Ma = zeros(1,imax);
    TT = zeros(1,imax);
    PT = zeros(1,imax);
    RT = zeros(1,imax);
    TS = zeros(1,imax);
    PS = zeros(1,imax);
    RS = zeros(1,imax);
    SS = zeros(1,imax);
    Ve = zeros(1,imax);

    Ma(1) = MIN;
    TT(1) = TTIN;
    PT(1) = PTIN;
    RT(1) = PT(1) / (RGAS * TT(1));

    TS(1) = TT(1) / (1 + 0.5*(gama-1)*Ma(1)^2);
    PS(1) = PT(1) / ((1 + 0.5*(gama-1)*Ma(1)^2)^(gama/(gama-1)));
    RS(1) = PS(1) / (RGAS * TS(1));
    SS(1) = sqrt(gama * RGAS * TS(1));
    Ve(1) = Ma(1) * SS(1);

    % march forward
    for i = 1:imax-1
        V = Ve(i);
        M = Ma(i);
        dA = A(i+1) - A(i);
        D = a_D(i);
        TTo = TT(i);

        if D <= 0
            error('Non-positive diameter at i=%d', i);
        end

        XM = 1 + 0.5*(gama-1)*M^2;
        XX = (M^2 - 1);

        % Avoid exact equality comparisons with floating point; use tolerance
        tol = 1e-8;
        if abs(XX) < tol
            % Mach near 1 -> singular behavior; stop with informative message
            warning('Mach ≈ 1 at i=%d (x=%.5f). Integration stopped to avoid singularity.', i, a_x(i));
            % Truncate arrays to current size and break
            Ma = Ma(1:i);
            Ve = Ve(1:i);
            TS = TS(1:i);
            PS = PS(1:i);
            RS = RS(1:i);
            SS = SS(1:i);
            TT = TT(1:i);
            PT = PT(1:i);
            a_x = a_x(1:i);
            a_D = a_D(1:i);
            break;
        else
            dV = (V/XX) * ((dA/A(i)) - XM * (dQ/CP/TTo) - 0.5*gama*M^2*(FR*dx/D));
        end

        Ve(i+1) = V + dV;

        Ro = RS(i);
        dRo = - Ro * (dV/V - dA/A(i));
        RS(i+1) = Ro + dRo;

        T = TS(i);
        dT = - T * ((gama-1)*(M^2)*(dV/V) - XM*(dQ/CP/TTo));
        TS(i+1) = T + dT;

        P = PS(i);
        dP = - P * (0.5*gama*M^2*(FR*dx/D) + gama*(M^2)*(dV/V));
        PS(i+1) = P + dP;

        SS(i+1) = sqrt(gama * RGAS * TS(i+1));
        Ma(i+1) = Ve(i+1) / SS(i+1);

        % recompute totals from local static
        TT(i+1) = TS(i+1) * (1 + 0.5*(gama-1)*Ma(i+1)^2);
        PT(i+1) = PS(i+1) * ((1 + 0.5*(gama-1)*Ma(i+1)^2)^(gama/(gama-1)));
        RT(i+1) = PT(i+1) / (RGAS * TT(i+1));
    end

    % store result
    results{icase} = struct('x',a_x,'D',a_D,'A',A,'Ma',Ma,'Ve',Ve,'TS',TS,'PS',PS,...
        'RS',RS,'SS',SS,'TT',TT,'PT',PT,'RT',RT,'FR',FR,'Q_total',Q_total);
end


% ---------------- PLOTTING ----------------
% Replace legends so they reflect actual cases
labels = {'Area change (isentropic)','Area + friction','Area + friction + heat'};

% Create output folder for results
output_folder = 'results';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Plot duct shape (diameter symmetric about centerline)
figure; hold on
plot(results{1}.x, results{1}.D, '-', 'LineWidth', 2.2);
plot(results{1}.x, -results{1}.D, '-', 'LineWidth', 2.2);
grid on;
xlabel('x (m)','fontsize',12); ylabel('Diameter (m)','fontsize',12);
title(sprintf('Convergent-divergent nozzle (A_{end}/A_{throat} ≈ %.2f)', ...
    results{1}.A(end)/min(results{1}.A)), 'fontsize',12);
legend('D(x)','Location','best');
hold off
saveas(gcf, fullfile(output_folder, 'nozzle_shape.png'));

% Helper for plotting and saving each variable
fn_plot = @(fieldname,ylabeltxt,figtitle,filename) ...
    plot_and_save(results, fieldname, labels, ylabeltxt, figtitle, output_folder, filename);

% Generate and save plots
fn_plot('Ve','Velocity (m/s)','Velocity Distribution','velocity.png');
fn_plot('TS','Static Temperature (K)','Static Temperature Distribution','static_temperature.png');
fn_plot('PS','Static Pressure (Pa)','Static Pressure Distribution','static_pressure.png');
fn_plot('RS','Density (kg/m^3)','Static Density Distribution','static_density.png');
fn_plot('SS','Speed of sound (m/s)','Speed of Sound Distribution','speed_of_sound.png');
fn_plot('Ma','Mach number','Mach Number Distribution','mach_number.png');
fn_plot('PT','Total Pressure (Pa)','Total Pressure Distribution','total_pressure.png');
fn_plot('TT','Total Temperature (K)','Total Temperature Distribution','total_temperature.png');
fn_plot('RT','Total density (kg/m^3)','Total Density Distribution','total_density.png');

% -------------------------------------------------------------------------
function plot_and_save(results, fieldname, labels, ylabeltxt, figtitle, folder, filename)
    figure; hold on; grid on;
    for k = 1:numel(results)
        r = results{k};
        if numel(r.x) == numel(r.(fieldname))
            plot(r.x, r.(fieldname), 'LineWidth', 2);
        else
            n = min([numel(r.x), numel(r.(fieldname))]);
            plot(r.x(1:n), r.(fieldname)(1:n), 'LineWidth', 2);
        end
    end
    xlabel('x (m)','fontsize',12);
    ylabel(ylabeltxt,'fontsize',12);
    title(figtitle,'fontsize',12);
    legend(labels,'Location','best');
    % Adjust axis padding for better visualization
    ax = gca;
    ax.XLim = [min(results{1}.x) max(results{1}.x)];
    ylim padded;
    hold off;
    saveas(gcf, fullfile(folder, filename));
end
