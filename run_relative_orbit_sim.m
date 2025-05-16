%% run_relative_orbit_sim.m
% ------------------------------------------------------------------------
% Relative-orbit tracking with Lyapunov feedback, 3D cartesian coordinates.
%
% Author : Jose Daniel Hoyos Giraldo (May-2025)
% ------------------------------------------------------------------------

clear; clc; close all;

% ========= USER-EDITABLE OPTIONS =======================================
pars = struct( ...
    ... % ---------- Chief (target) orbit -------------------------------
    "a_c",      15000 , ...   % km  (semi-major axis)
    "e_c",           0.3 , ... % –   (eccentricity)
    "i_c",            30 , ... % deg (inclination)
    "RAAN_c",          0 , ... % deg (right ascension of ascending node)
    "omega_c",         0 , ... % deg (argument of perigee)
    "nu0_c",          45 , ... % deg (true anomaly at t = 0)
    ... % ---------- Deputy initial (natural) orbit ---------------------
    "a_d",      11000 , ...   % km
    "e_d",           0.3 , ... % –
    "i_d",            25 , ... % deg
    "RAAN_d",          0 , ... % deg
    "omega_d",         0 , ... % deg
    "nu0_d",          45 , ... % deg
    ... % ---------- Control & simulation settings ----------------------
    "KrScale",        2 , ... % (choose 0.4 | 2 | 10 …) for underdamped, critically damped, overdamped
    "u_max",      1000 , ...  % km/s²   (physical thrust saturation) Let a high value for no max constraint
    "tf_hours",       8 , ... % total simulation time (hours)
    "dt_sec",        0.1 , ... % fixed ODE step (seconds, real time)
    ... % ---------- Visualisation flags -------------------------------
    "animate",      true , ... % true  → dual-panel animation
    "fullPlots",    true   ... % true  → time-history figures (already inside engine)
);
% =======================================================================

% ---- RUN the simulation ------------------------------------------------
out = simulateRelativeOrbit(pars);

% ---- Quick console summary --------------------------------------------
fprintf('\n=== Simulation complete ===\n');
fprintf('Simulated %g hours with γ = %g, thrust sat = %.3g km/s²\n', ...
        pars.tf_hours, pars.KrScale, pars.u_max);
fprintf('Total control effort  ∫||u|| dt  = %.4f  km/s\n', out.J_u);