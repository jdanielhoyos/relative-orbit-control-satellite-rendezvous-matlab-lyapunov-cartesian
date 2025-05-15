function out = simulateRelativeOrbit(varargin)
% USAGE 1:  simulateRelativeOrbit(optsStruct)
% USAGE 2:  simulateRelativeOrbit("field1",value1,"field2",value2,…)
%
%            recognised fields (defaults in brackets)
%   a_c,e_c,i_c,RAAN_c,omega_c,nu0_c     Chief OE  [10000,0.3,30,0,0,45]
%   a_d,e_d,i_d,RAAN_d,omega_d,nu0_d     Deputy OE [11000,0.3,25,0,0,45]
%   KrScale      – γ in P = γ·diag(√Kr)              [20]
%   u_max        – thrust sat (km/s², physical)      [0.001]
%   tf_hours     – simulation length (hours)         [10]
%   dt_sec       – fixed ODE step (seconds)          [0.1]
%   animate      – true/false; run the big animation [true]
% --------------------------------------------------------------
% Outputs (structure “out”) – time axis, histories, effort, etc.
% --------------------------------------------------------------
% Author : Jose Daniel Hoyos Giraldo  (May-2025)
% --------------------------------------------------------------

%% -------- 0.  Option parser -------------------------------------------
def = struct( ...
  "a_c",10000,"e_c",0.3,"i_c",30,"RAAN_c",0,"omega_c",0,"nu0_c",45, ...
  "a_d",11000,"e_d",0.3,"i_d",25,"RAAN_d",0,"omega_d",0,"nu0_d",45, ...
  "KrScale",20,"u_max",0.001,"tf_hours",10,"dt_sec",0.1,"animate",true );

if nargin==1 && isstruct(varargin{1})
    user = varargin{1};
    f = fieldnames(user);  for k=1:numel(f),  def.(f{k}) = user.(f{k}); end
else
    if mod(nargin,2)~=0,  error("Name/value pairs required.");  end
    for k=1:2:nargin
        name = varargin{k};  val = varargin{k+1};
        if ~isfield(def,name), error("Unknown option '%s'",name); end
        def.(name) = val;
    end
end
p = def;                        % shorter alias

%% -------- 1.  Constants / initial states ------------------------------
mu = 398600.4418;               % km^3/s^2
L  = p.a_c;                     % characteristic length (chief SMA)
T  = sqrt(L^3/mu);              % characteristic time

% Absolute states at t=0
[rc0, vc0] = orbitalElements2State(p.a_c,p.e_c,p.i_c,p.RAAN_c,p.omega_c,p.nu0_c,mu);
[rd0, vd0] = orbitalElements2State(p.a_d,p.e_d,p.i_d,p.RAAN_d,p.omega_d,p.nu0_d,mu);

x0_dim = [rd0 - rc0;  vd0 - vc0];          % relative (km, km/s)
x0     = [x0_dim(1:3)/L;  x0_dim(4:6)/(L/T)];  % normalised

%% -------- 2.  Control gains & limits ----------------------------------
Kr = eye(3);
P  = diag(p.KrScale*sqrt(diag(Kr)));
u_max_n = p.u_max / (L/T^2);               % normalised

%% -------- 3.  Integration ---------------------------------------------
Tf_n    = p.tf_hours*3600 / T;
dt_n    = p.dt_sec / T;
tspan_n = 0:dt_n:Tf_n;

f = @(t,x) nonlinear_relative_dynamics(t,x,Kr,P,u_max_n, ...
         p.a_c,p.e_c,p.i_c,p.RAAN_c,p.omega_c,deg2rad(p.nu0_c), ...
         mu,L,T);

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t_out,x_out] = ode45(f,tspan_n,x0,opts);

delta_r_n = x_out(:,1:3);  delta_v_n = x_out(:,4:6);

%% -------- 4.  Control history (post-process, identical to your code) ---
mu_n = 1;   u_hist_n = zeros(size(delta_r_n));
for k = 1:length(t_out)
    dr = delta_r_n(k,:)';  dv = delta_v_n(k,:)';
    t_phys_k = t_out(k)*T;
    [~,~,nu_ck] = propagateKeplerOrbit(p.a_c,p.e_c,deg2rad(p.nu0_c),mu,t_phys_k);
    [r_ck,~] = orbitalElements2State(p.a_c,p.e_c,p.i_c,p.RAAN_c,...
                                     p.omega_c,rad2deg(nu_ck),mu);
    r_c_norm_k = r_ck/L;
    r_deputy_k = dr + r_c_norm_k;
    delta_a_k  = -mu_n*( r_deputy_k/norm(r_deputy_k)^3 - r_c_norm_k/norm(r_c_norm_k)^3 );
    u_k = -Kr*dr - delta_a_k - P*dv;
    if norm(u_k)>u_max_n, u_k = u_max_n*u_k/norm(u_k); end
    u_hist_n(k,:) = u_k';
end

%% -------- 5.  Denormalise for plotting (km, km/s, hr) -----------------
time_hr = t_out*T/3600;
delta_r = delta_r_n*L;
delta_v = delta_v_n*(L/T);
u_hist  = u_hist_n*(L/T^2);

%% -------- 6.  ORIGINAL PLOTS (unchanged code begins) ------------------
% *** EVERYTHING below is copied verbatim from your script ***

%% Relative Position
figure('Color','w');
plot(time_hr,delta_r); hold on;
plot(time_hr,vecnorm(delta_r,2,2),'k--','LineWidth',1.5);
xlabel('Time [hr]'); ylabel('\delta r [km]');
title('Relative Position Error');
legend({'\delta r_1','\delta r_2','\delta r_3','||\delta r||'},'Location','best');
grid on;

%% Relative Velocity
figure('Color','w');
plot(time_hr,delta_v); hold on;
plot(time_hr,vecnorm(delta_v,2,2),'k--','LineWidth',1.5);
xlabel('Time [hr]'); ylabel('\delta v [km/s]');
title('Relative Velocity Error');
legend({'\delta v_1','\delta v_2','\delta v_3','||\delta v||'},'Location','best');
grid on;

%% Control Input
figure('Color','w');
plot(time_hr,u_hist); hold on;
plot(time_hr,vecnorm(u_hist,2,2),'k--','LineWidth',1.5);
xlabel('Time [hr]'); ylabel('u [km/s^2]');
title('Control Input');
legend({'u_1','u_2','u_3','||u||'},'Location','best');
grid on;

%% Combined Error (log scale)
Delta_combined = sqrt((delta_r/L).^2 + (delta_v*T/L).^2);
Delta_mag = vecnorm(Delta_combined,2,2);
figure('Color','w');
semilogy(time_hr,Delta_combined); hold on;
semilogy(time_hr,Delta_mag,'k--','LineWidth',1.5);
xlabel('Time [hr]'); ylabel('\Delta [km] (log scale)');
title('Combined Position-Velocity Error');
legend({'\Delta_1','\Delta_2','\Delta_3','||\Delta||'},'Location','best');
grid on;

%% 3-D Relative Trajectory
figure;
plot3(delta_r(:,1),delta_r(:,2),delta_r(:,3),'b','DisplayName','Path'); hold on;
plot3(delta_r(1,1),delta_r(1,2),delta_r(1,3),'go','MarkerSize',10,'DisplayName','Initial');
plot3(delta_r(end,1),delta_r(end,2),delta_r(end,3),'ro','MarkerSize',10,'DisplayName','Final');
plot3(0,0,0,'kx','MarkerSize',12,'LineWidth',2,'DisplayName','Target');
xlabel('\delta r_1 [km]'); ylabel('\delta r_2 [km]'); zlabel('\delta r_3 [km]');
title('3D Relative Trajectory'); grid on; axis equal; axis square; legend;
set(gcf,'Color','k'); set(gca,'Color','k');
set(gca,'XColor','w','YColor','w','ZColor','w');
set(gca,'GridColor','w','MinorGridColor','w');
set(get(gca,'Title'),'Color','w');
set(legend,'TextColor','w','Color','k');
lines = findobj(gcf,'Type','Line');
set(lines,'Color','w'); set(lines(2),'MarkerEdgeColor','y');

%% Absolute Orbits (Comparison)  --- identical block --------------------
chief_abs  = zeros(length(t_out),3);
deputy_abs = zeros(length(t_out),3);
chief_abs(1,:)  = rc0';   deputy_abs(1,:) = rd0';

for k = 2:length(t_out)
    t_phys_k = t_out(k)*T;
    [~,~,nu_ck] = propagateKeplerOrbit(p.a_c,p.e_c,deg2rad(p.nu0_c),mu,t_phys_k);
    [r_ck,~] = orbitalElements2State(p.a_c,p.e_c,p.i_c,p.RAAN_c,p.omega_c,rad2deg(nu_ck),mu);
    chief_abs(k,:) = r_ck';
    deputy_abs(k,:) = r_ck' + delta_r(k,:);
end

theta_ref = linspace(0,2*pi,500);
initial_deputy_orbit = zeros(3,numel(theta_ref));
target_orbit          = zeros(3,numel(theta_ref));
for j = 1:numel(theta_ref)
    [r_id,~] = orbitalElements2State(p.a_d,p.e_d,p.i_d,p.RAAN_d,p.omega_d,rad2deg(theta_ref(j)),mu);
    [r_tg,~] = orbitalElements2State(p.a_c,p.e_c,p.i_c,p.RAAN_c,p.omega_c,rad2deg(theta_ref(j)),mu);
    initial_deputy_orbit(:,j) = r_id;  target_orbit(:,j) = r_tg;
end
figure;
plot_Earth(500,4); hold on;
plot3(chief_abs(:,1),chief_abs(:,2),chief_abs(:,3),'-','Color',[0.6 0.6 1], ...
      'LineWidth',1.2,'DisplayName','Target Orbit (Chief)');
plot3(initial_deputy_orbit(1,:),initial_deputy_orbit(2,:),initial_deputy_orbit(3,:), ...
      '-','Color',[0.4 1 0.4],'LineWidth',1.2,'DisplayName','Initial Deputy Orbit');
plot3(deputy_abs(:,1),deputy_abs(:,2),deputy_abs(:,3),'-','Color','white', ...
      'LineWidth',1.2,'DisplayName','Deputy (Controlled)');
plot3(chief_abs(1,1),chief_abs(1,2),chief_abs(1,3),'oc','MarkerSize',8,'MarkerFaceColor','c',...
      'DisplayName','Chief Start');
plot3(chief_abs(end,1),chief_abs(end,2),chief_abs(end,3),'+b','MarkerSize',15,'LineWidth',3,...
      'DisplayName','Chief Final');
plot3(deputy_abs(1,1),deputy_abs(1,2),deputy_abs(1,3),'oy','MarkerSize',8,'MarkerFaceColor','y',...
      'DisplayName','Deputy Start');
plot3(deputy_abs(end,1),deputy_abs(end,2),deputy_abs(end,3),'xr','MarkerSize',15,'LineWidth',3,...
      'DisplayName','Deputy Final');
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('Orbits and Controlled Deputy Trajectory');
grid on; axis equal; legend('Location','best');
set(gcf,'Color','k'); set(gca,'Color','k');
set(gca,'XColor','w','YColor','w','ZColor','w','GridColor','w','MinorGridColor','w');
set(get(gca,'Title'),'Color','w'); set(legend,'TextColor','w','Color','k');

%% Control effort
time_sec = time_hr*3600;  u_norm = vecnorm(u_hist,2,2);
integral_u = trapz(time_sec,u_norm);
fprintf('Total control effort (∫||u|| dt) = %.4f km/s\n',integral_u);

%% Big animation (identical) -------------------------------------------
if p.animate
figure('WindowState','maximized','Color','k','MenuBar','none','ToolBar','none');
trail_length = 100;
ax_abs = subplot(1,2,1); hold(ax_abs,'on'); axis(ax_abs,'equal'); axis(ax_abs,'off');
title(ax_abs,'Absolute Orbits','Color','w');
ax_rel = subplot(1,2,2); hold(ax_rel,'on'); axis(ax_rel,'equal'); grid(ax_rel,'on');
xlabel(ax_rel,'\delta r_1 [km]','Color','w');
ylabel(ax_rel,'\delta r_2 [km]','Color','w');
zlabel(ax_rel,'\delta r_3 [km]','Color','w');
title(ax_rel,'Relative Trajectory','Color','w');
step       = 3000;
time_anim  = time_hr(1:step:end);
chief_anim = chief_abs(1:step:end,:);
deputy_anim= deputy_abs(1:step:end,:);
rel_anim   = delta_r(1:step:end,:);
axes(ax_abs);
view(40.9514,16.3231)
gh = plot_Earth(500,4);
h_chief  = scatter3(chief_anim(1,1),chief_anim(1,2),chief_anim(1,3),'x','MarkerEdgeColor',[1 0 0],...
                    'SizeData',80,'LineWidth',2);
h_deputy = scatter3(deputy_anim(1,1),deputy_anim(1,2),deputy_anim(1,3),'+','MarkerEdgeColor',[1 1 1],...
                    'SizeData',80,'LineWidth',2);
h_time_abs = text(0.02,0.95,'Time = 0.0 hr','Units','normalized','Color','w','FontSize',14);
plot3(ax_abs,target_orbit(1,:),target_orbit(2,:),target_orbit(3,:),'--r');
plot3(ax_abs,initial_deputy_orbit(1,:),initial_deputy_orbit(2,:),initial_deputy_orbit(3,:),...
      '--','Color',[0.4 1 0.4]);
h_deputy_trail = plot3(deputy_anim(1,1),deputy_anim(1,2),deputy_anim(1,3),'-','Color',[1 1 1], ...
                       'LineWidth',1,'HandleVisibility','off');
axes(ax_rel);
plot3(ax_rel,rel_anim(:,1),rel_anim(:,2),rel_anim(:,3),'--','LineWidth',1);
scatter3(ax_rel,rel_anim(1,1),rel_anim(1,2),rel_anim(1,3),80,'o',...
         'MarkerEdgeColor','g','MarkerFaceColor','g');
scatter3(ax_rel,0,0,0,100,'x','LineWidth',2,'MarkerEdgeColor','r');
h_rel       = scatter3(ax_rel,rel_anim(1,1),rel_anim(1,2),rel_anim(1,3),60,'filled',...
                       'MarkerEdgeColor','w','MarkerFaceColor','w');
h_rel_trail = plot3(ax_rel,rel_anim(1,1),rel_anim(1,2),rel_anim(1,3),'-','Color','w','LineWidth',1);
lg_rel = legend(ax_rel,{'Trajectory','Initial','Target'},'Location','best');
set(lg_rel,'TextColor','w','Color','k','Box','on');
for ax = [ax_abs, ax_rel]
    set(ax,'Color','k','XColor','w','YColor','w','ZColor','w','GridColor','w','MinorGridColor','w');
end
trail_len = 10; pause_dur = 0.005;
for k = 1:length(time_anim)
    if k==1, pause(1), end
    axes(ax_abs);
    
    delete(gh);
    gh = plot_Earth(50,mod(time_anim(k)*15,360));
    set(h_chief,'XData',chief_anim(k,1),'YData',chief_anim(k,2),'ZData',chief_anim(k,3));
    set(h_deputy,'XData',deputy_anim(k,1),'YData',deputy_anim(k,2),'ZData',deputy_anim(k,3));
    set(h_time_abs,'String',sprintf('Time = %.1f hr',time_anim(k)));
    % start_idx = max(1,k-trail_length); cur = start_idx:k; %If want trail with only last trail_len
    cur = 1:k;
    set(h_deputy_trail,'XData',deputy_anim(cur,1),'YData',deputy_anim(cur,2),'ZData',deputy_anim(cur,3));
    axes(ax_rel);
    set(h_rel,'XData',rel_anim(k,1),'YData',rel_anim(k,2),'ZData',rel_anim(k,3));
    %idx0 = max(1,k-trail_len); %If want trail with only last trail_len
    %steps
    idx0 = 1;
    set(h_rel_trail,'XData',rel_anim(idx0:k,1),'YData',rel_anim(idx0:k,2),'ZData',rel_anim(idx0:k,3));
    drawnow limitrate; pause(pause_dur);
end
end % animate

%% -------- 7.  Return structure ----------------------------------------
out.time_hr    = time_hr;
out.delta_r    = delta_r;
out.delta_v    = delta_v;
out.u_hist     = u_hist;
out.J_u        = integral_u;
out.chief_abs  = chief_abs;
out.deputy_abs = deputy_abs;

% ---- END MAIN FUNCTION -------------------------------------------------
end



% ======================  LOCAL SUB-FUNCTIONS  ==========================

function dx = relDyn(t,x,Kr,P,u_max_n,p,L,T,mu)
    [delta_a, r_c_norm] = inertialTerms(t,x(1:3),p,L,T,mu);

    % Control law with saturation
    u = -Kr*x(1:3) - delta_a - P*x(4:6);
    if norm(u) > u_max_n,  u = u_max_n*u/norm(u);  end

    dx = [x(4:6);
          delta_a + u];
end
% -----------------------------------------------------------------------
function [delta_a, r_c_norm] = inertialTerms(t,delta_r_n,p,L,T,mu)
    % Chief position at time t (normalised units inside)
    t_sec = t*T;
    [~,~,nu_c] = propagateKeplerOrbit(...
         p.a_c, p.e_c, deg2rad(p.nu0_c), mu, t_sec);
    [r_c_phys,~] = orbitalElements2State(...
         p.a_c, p.e_c, p.i_c, p.RAAN_c, p.omega_c, rad2deg(nu_c), mu);
    r_c_norm = r_c_phys./L;

    r_deputy = delta_r_n(:) + r_c_norm;
    mu_n = 1;
    delta_a = -mu_n*( r_deputy/norm(r_deputy)^3 - r_c_norm/norm(r_c_norm)^3 );
end
% -----------------------------------------------------------------------
function [chief_abs,deputy_abs] = buildAbsoluteTraj(t_out,T,p,delta_r,mu)
    N = numel(t_out);
    chief_abs  = zeros(N,3);
    deputy_abs = zeros(N,3);

    for k = 1:N
        t_sec = t_out(k)*T;
        [~,~,nu_ck] = propagateKeplerOrbit(p.a_c,p.e_c,deg2rad(p.nu0_c),mu,t_sec);
        [r_ck,~] = orbitalElements2State(p.a_c,p.e_c,p.i_c,p.RAAN_c,...
                                         p.omega_c,rad2deg(nu_ck),mu);
        chief_abs(k,:) = r_ck.';

        deputy_abs(k,:) = r_ck.' + delta_r(k,:);
    end
end
% -----------------------------------------------------------------------
function plotTimeHistories(t,dr,dv,u,L,T)
    figure('Color','w');
    subplot(3,1,1);
    plot(t,dr); hold on; plot(t,vecnorm(dr,2,2),'k--','LineWidth',1.4);
    xlabel('Time [h]'); ylabel('\delta r [km]'); title('Relative position');

    subplot(3,1,2);
    plot(t,dv); hold on; plot(t,vecnorm(dv,2,2),'k--','LineWidth',1.4);
    xlabel('Time [h]'); ylabel('\delta v [km/s]'); title('Relative velocity');

    subplot(3,1,3);
    plot(t,u);  hold on; plot(t,vecnorm(u,2,2),'k--','LineWidth',1.4);
    xlabel('Time [h]'); ylabel('u [km/s^2]'); title('Control input');
end
% -----------------------------------------------------------------------
function animateTraj(time_hr,chief,deputy,dr,p)
    % Dual-panel animation that re-uses the user's plot_Earth() utility.
    %
    % Requires a file  plot_Earth(radius_km, GST_deg)  in the MATLAB path.
    % The first argument controls the sphere size used for plotting
    % (you were using 500 km and 50 km in your original script).
    % The second argument is the Greenwich sidereal angle that simply
    % spins the texture so the Earth appears to rotate.

    figure('Color','k','WindowState','maximized');
    step  = 3000;               % draw every Nth sample
    axAbs = subplot(1,2,1); hold on; axis equal off;
    axRel = subplot(1,2,2); hold on; axis equal; grid on;
    title(axAbs,'Absolute orbits','Color','w');
    title(axRel,'Relative trajectory','Color','w');
    set(axRel,'Color','k','XColor','w','YColor','w','ZColor','w');

    % ----- initial draw -------------------------------------------------
    gh = plot_Earth(50,0);      % <--- YOUR OWN FUNCTION, unchanged
    plot3(axAbs,chief(:,1),chief(:,2),chief(:,3),'-','Color',[.6 .6 1]);
    plot3(axRel,dr(:,1),dr(:,2),dr(:,3),'--','Color',[.6 .6 .6]);

    hC = plot3(axAbs,chief(1,1),chief(1,2),chief(1,3),'xr','LineWidth',2);
    hD = plot3(axAbs,deputy(1,1),deputy(1,2),deputy(1,3),'+w','LineWidth',2);
    hR = plot3(axRel,dr(1,1),   dr(1,2),   dr(1,3),  'ow','MarkerFaceColor','w');

    % ----- animation loop ----------------------------------------------
    for k = 1:step:numel(time_hr)
        delete(gh);                       % spin the globe
        gh = plot_Earth(50,mod(time_hr(k)*15,360));

        set(hC,'XData',chief(k,1),'YData',chief(k,2),'ZData',chief(k,3));
        set(hD,'XData',deputy(k,1),'YData',deputy(k,2),'ZData',deputy(k,3));
        set(hR,'XData',dr(k,1),   'YData',dr(k,2),   'ZData',dr(k,3));

        drawnow limitrate;
    end
end
% -----------------------------------------------------------------------
function h = plotEarth(ax,dalpha,GSTdeg)
    % Very lightweight textured Earth (merely a sphere w/ basic lighting)
    R_E = 6378;  [X,Y,Z] = sphere(120);
    h = surf(ax,R_E*X,R_E*Y,R_E*Z, 'EdgeColor','none', ...
        'FaceColor',[0.05 0.2 0.6], 'FaceLighting','gouraud');
    camlight(ax,'headlight'); material(ax,'dull');
    rotate(h,[0 0 1],GSTdeg);   % spin by Greenwich sidereal angle
    alpha(h,1-dalpha/100);      % simple daylight‐night terminator
end
% -----------------------------------------------------------------------
function [M,E,nu] = propagateKeplerOrbit(a,e,nu0,mu,t)
    E0 = 2*atan( sqrt((1-e)/(1+e))*tan(nu0/2) );
    M0 = E0 - e*sin(E0);  n = sqrt(mu/a^3);
    M  = mod(M0 + n*t, 2*pi);

    % Newton–Raphson for Kepler’s equation
    E = M;  for i=1:50
        dE = (E - e*sin(E) - M)/(1 - e*cos(E));
        E  = E - dE;  if abs(dE) < 1e-12, break; end
    end
    nu = 2*atan( sqrt((1+e)/(1-e))*tan(E/2) );
end
% -----------------------------------------------------------------------
function [r,v] = orbitalElements2State(a,e,i,RAAN,omega,nu,mu)
    % Angle conversion
    i=deg2rad(i); RAAN=deg2rad(RAAN); omega=deg2rad(omega); nu=deg2rad(nu);

    p = a*(1-e^2);
    r_pf = p/(1+e*cos(nu))*[cos(nu); sin(nu); 0];
    v_pf = sqrt(mu/p)*[-sin(nu); e+cos(nu); 0];

    R3_W = [cos(RAAN),-sin(RAAN),0; sin(RAAN),cos(RAAN),0; 0 0 1];
    R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_w = [cos(omega),-sin(omega),0; sin(omega),cos(omega),0; 0 0 1];

    Q = R3_W*R1_i*R3_w;
    r = Q*r_pf;  v = Q*v_pf;
end
% =======================================================================
function dx = nonlinear_relative_dynamics(t, x, Kr, P, u_max, ...
                                          a_c, e_c, i_c, RAAN_c, ...
                                          omega_c, nu_c, mu, L, T)
    % Convert normalized time -> real seconds
    t_phys = t * T;

    % Chief's true anomaly at this time
    [~, ~, nu_chief] = propagateKeplerOrbit(a_c, e_c, nu_c, mu, t_phys);

    % Chief's current position in ECI (phys. units), then normalize
    [r_c_phys, ~] = orbitalElements2State(a_c, e_c, i_c, RAAN_c, ...
                                          omega_c, rad2deg(nu_chief), mu);
    r_c_norm = r_c_phys / L;

    % Extract deputy's relative state
    delta_r = x(1:3);
    delta_v = x(4:6);

    % Differential gravitational acceleration
    r_deputy = delta_r + r_c_norm;
    mu_n = 1;  % normalized gravitational parameter
    delta_a = -mu_n * ( r_deputy / norm(r_deputy)^3 ...
                        - r_c_norm / norm(r_c_norm)^3 );

    % Control law
    u = -Kr*delta_r - delta_a - P*delta_v;

    % Saturation
    if norm(u) > u_max
        u = (u_max / norm(u)) * u;
    end

    % State derivative (normalized)
    dx = zeros(6,1);
    dx(1:3) = delta_v;
    dx(4:6) = delta_a + u;
end