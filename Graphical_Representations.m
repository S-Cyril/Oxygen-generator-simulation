function Graphical_Representations

R      = 8.314;    T_in = 298.0;   Pbar = 4.5;  P_in = Pbar*1e5;
yN2_in = 0.79;

L     = 0.55;   eps = 0.37;   rho_p = 1100.0;   dp = 2e-3;
u_in  = 0.09;   Dax = 7.5e-5;  kLDF = 0.11;     mu = 1.8e-5;
lambda_eff = 0.15;  Cp_g_m = 29.0;  Cp_s = 900.0;
dH_N2 = 2.0e4;   dH_O2 = 1.2e4;

qmaxN2 = 1.8;      bN2bar = 0.07;
qmaxO2 = 0.45;     bO2bar = 0.02;
bN2 = bN2bar/1e5;  bO2 = bO2bar/1e5;

Nx  = 80;  x  = linspace(0,L,Nx).';  dx = L/(Nx-1);
tend = 60.0;   

yN2_0 = zeros(Nx,1);   qN2_0 = zeros(Nx,1);
qO2_0 = zeros(Nx,1);   T0    = T_in*ones(Nx,1);
y0    = [yN2_0; qN2_0; qO2_0; T0];


pars = struct('Nx',Nx,'dx',dx,'eps',eps,'rho_p',rho_p,'R',R, ...
              'Dax',Dax,'kLDF',kLDF,'qmaxN2',qmaxN2,'qmaxO2',qmaxO2, ...
              'bN2',bN2,'bO2',bO2,'yN2_in',yN2_in, ...
              'dp',dp,'mu',mu,'Cp_g_m',Cp_g_m,'Cp_s',Cp_s, ...
              'lambda_eff',lambda_eff,'dH_N2',dH_N2,'dH_O2',dH_O2, ...
              'P_in',P_in,'u_in',u_in,'T_in',T_in);


opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.25);
[Tsol, Y] = ode15s(@(t,y) rhs_yN2_xT_fw(t,y,pars), [0 tend], y0, opts);


yN2_out_trace = Y(:, 1:pars.Nx);
yN2_out_trace = yN2_out_trace(:, end);
thr10 = 0.10*yN2_in;  thr50 = 0.50*yN2_in;

tB10  = localFirstCross_fw(yN2_out_trace, Tsol, thr10);
tB50  = localFirstCross_fw(yN2_out_trace, Tsol, thr50);


snapT = [0, 10, 20, 40, tend];   
idx   = arrayfun(@(ts) find(abs(Tsol-ts)==min(abs(Tsol-ts)),1,'first'), snapT);


yN2_end = Y(end,1:Nx).';
T_end   = Y(end,3*Nx+1:4*Nx).';
prof_end = compute_pressure_velocity_x_fw(pars, yN2_end, T_end);
P_end    = prof_end.P;     
Ctot_end = prof_end.Ctot;  


figure('Name','Temperature profile'); hold on; grid on;
for m = 1:numel(idx)
    T_snap = Y(idx(m), 3*Nx+1:4*Nx).';
    plot(x, T_snap, 'LineWidth', 1.5);
end
xlabel('Axial position, x (m)');
ylabel('Temperature, T (K)');
legend(arrayfun(@(ts) sprintf('t = %.0f s', ts), snapT, 'UniformOutput', false), ...
       'Location','best');
title('Temperature profile along the bed');


figure('Name','Pressure drop'); grid on;
plot(x, P_end/1e5, 'LineWidth', 1.8);
xlabel('Axial position, x (m)');
ylabel('Pressure, P (bar)');
dPbar = (P_end(1)-P_end(end))/1e5;
title(sprintf('Pressure profile at final time (\\DeltaP = %.2f bar)', dPbar));


figure('Name','Species profile'); hold on; grid on;
yN2_f = yN2_end;
yO2_f = 1 - yN2_f;
plot(x, yN2_f, 'LineWidth', 1.8);
plot(x, yO2_f, 'LineWidth', 1.8);
xlabel('Axial position, x (m)');
ylabel('Mole fraction, y');
legend('y_{N_2}(x)','y_{O_2}(x)','Location','best');
title('Composition from inlet (x=0) to outlet (x=L) at final time');


figure('Name','Breakthrough'); hold on; grid on;
plot(Tsol, yN2_out_trace, 'LineWidth', 1.8);
yl = ylim;
if ~isnan(tB10)
    plot([tB10 tB10], yl, '--', 'LineWidth', 1.0);
end
if ~isnan(tB50)
    plot([tB50 tB50], yl, '--', 'LineWidth', 1.0);
end
xlabel('Time, t (s)'); ylabel('Outlet y_{N_2} (–)');
leg = {'y_{N_2,out}(t)'};
if ~isnan(tB10), leg{end+1}=sprintf('t_{B10}=%.2f s', tB10); end
if ~isnan(tB50), leg{end+1}=sprintf('t_{B50}=%.2f s', tB50); end
legend(leg,'Location','best');
title('Breakthrough behaviour at x = L (air in at x = 0)');

end



function dy = rhs_yN2_xT_fw(~,y,par)
    Nx = par.Nx; dx = par.dx; eps = par.eps; rho_p = par.rho_p;
    kLDF = par.kLDF;
    qmaxN2 = par.qmaxN2; qmaxO2 = par.qmaxO2; bN2 = par.bN2; bO2 = par.bO2;
    Dax = par.Dax; Cp_g_m = par.Cp_g_m; Cp_s = par.Cp_s; lambda_eff = par.lambda_eff;
    dH_N2 = par.dH_N2; dH_O2 = par.dH_O2;  yN2_in = par.yN2_in;

    
    yN2 = y(1:Nx);
    qN2 = y(Nx+1:2*Nx);
    qO2 = y(2*Nx+1:3*Nx);
    T   = y(3*Nx+1:4*Nx);

   
    prof = compute_pressure_velocity_x_fw(par, yN2, T);
    P    = prof.P;    U = prof.U;    Ctot = prof.Ctot;

    
    yN2 = max(1e-9, min(1-1e-9, yN2));   yO2 = 1 - yN2;

    
    PN2 = yN2 .* P;     PO2 = yO2 .* P;
    den = 1 + bN2.*PN2 + bO2.*PO2;
    qN2s = qmaxN2 * (bN2.*PN2) ./ den;
    qO2s = qmaxO2 * (bO2.*PO2) ./ den;

    dqN2dt = kLDF * (qN2s - qN2);
    dqO2dt = kLDF * (qO2s - qO2);

    
    d2y = zeros(Nx,1);
    d2y(2:end-1) = (yN2(3:end) - 2*yN2(2:end-1) + yN2(1:end-2))/dx^2;
    d2y(1)   = (yN2(2) - yN2(1))/dx^2;
    d2y(end) = (yN2(end-1) - yN2(end))/dx^2;

    dy_dx = zeros(Nx,1);
    dy_dx(1) = (yN2(1) - yN2_in)/dx;    
    dy_dx(2:end) = (yN2(2:end) - yN2(1:end-1))/dx;

    dyN2dt = ( Dax.*Ctot .* d2y - U.*Ctot .* dy_dx - ((1-eps)/eps)*rho_p.*dqN2dt ) ./ (eps.*Ctot);

   
    Ce   = eps.*Ctot.*Cp_g_m + (1-eps).*rho_p.*Cp_s;

    d2T = zeros(Nx,1);
    d2T(2:end-1) = (T(3:end) - 2*T(2:end-1) + T(1:end-2))/dx^2;
    d2T(1)   = (T(2) - T(1))/dx^2;
    d2T(end) = (T(end-1) - T(end))/dx^2;

    dT_dx = zeros(Nx,1);
    dT_dx(1) = 0;
    dT_dx(2:end) = (T(2:end) - T(1:end-1))/dx;

    Qads = (1-eps).*rho_p .* ( dH_N2.*dqN2dt + dH_O2.*dqO2dt );

    dTdt = ( lambda_eff.*d2T + Qads - (U.*Ctot*Cp_g_m).*dT_dx ) ./ Ce;

    dy = [dyN2dt; dqN2dt; dqO2dt; dTdt];
end

function prof = compute_pressure_velocity_x_fw(par, yN2, T)
    Nx = par.Nx; dx = par.dx; R = par.R; eps = par.eps;
    dp = par.dp; mu = par.mu; P_in = par.P_in; u_in = par.u_in;

    yN2 = max(1e-9, min(1-1e-9, yN2));  yO2 = 1 - yN2;
    M_N2 = 0.0280134; M_O2 = 0.031998;

   
    Ctot_in = P_in/(R*T(1));
    Gm = u_in * Ctot_in;  

    P   = zeros(Nx,1);  U   = zeros(Nx,1);  Ctt = zeros(Nx,1);
    P(1) = P_in;

    for i = 1:Nx-1
        Ctt(i) = P(i)/(R*T(i));
        U(i)   = Gm / Ctt(i);
        rho_g  = Ctt(i) * (yN2(i)*M_N2 + (1-yN2(i))*M_O2);

        visc  = 150*((1-eps)^2/eps^3) * mu * U(i) / (dp^2);
        inert = 1.75*((1-eps)/eps^3) * rho_g * U(i)^2 / dp;
        dPdx  = -(visc + inert);

        P(i+1) = max(1e3, P(i) + dPdx*dx); 
    end
    Ctt(end) = P(end)/(R*T(end));
    U(end)   = Gm / Ctt(end);

    prof = struct('P',P,'U',U,'Ctot',Ctt);
end

function tB = localFirstCross_fw(ytrace, tvec, thr)
    idx = find(ytrace >= thr, 1, 'first');
    if isempty(idx)
        tB = NaN;
    elseif idx == 1
        tB = tvec(1);
    else
        y0c = ytrace(idx-1); y1c = ytrace(idx);
        t0c = tvec(idx-1);   t1c = tvec(idx);
        if y1c == y0c
            tB = t1c;
        else
            tB = t0c + (thr - y0c)*(t1c - t0c)/(y1c - y0c);
        end
    end
end
