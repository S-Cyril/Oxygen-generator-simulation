

function Oxygen_Purity_Test
    
    R     = 8.314;              
    T_in  = 298.0;              
    Pbar  = 4.5;                
    P_in  = Pbar*1e5;           
    yN2_in = 0.79;              

    
    L    = 0.55;                
    eps  = 0.37;                
    rho_p = 1100.0;             
    dp   = 2e-3;                
    
    u_in = 0.09;                
    Dax  = 7.5e-5;              
    kLDF = 0.11;                
    mu   = 1.8e-5;              
    lambda_eff = 0.15;          
    Cp_g_m = 29.0;              
    Cp_s   = 900.0;             
    dH_N2  = 2.0e4;             
    dH_O2  = 1.2e4;             

    
    qmaxN2 = 1.8;               
    bN2bar = 0.07;              
    qmaxO2 = 0.45;              
    bO2bar = 0.02;              
    bN2 = bN2bar/1e5;           
    bO2 = bO2bar/1e5;           

    
    Nx   = 80;
    x    = linspace(0,L,Nx).'; 
    dx   = L/(Nx-1);
    tend = 18.0;                

   
    yN2_0 = yN2_in*ones(Nx,1);
    qN2_0 = zeros(Nx,1);
    qO2_0 = zeros(Nx,1);
    T0    = T_in*ones(Nx,1);
    y0    = [yN2_0; qN2_0; qO2_0; T0];

    
    pars = struct('Nx',Nx,'dx',dx,'eps',eps,'rho_p',rho_p,'R',R, ...
                  'Dax',Dax,'kLDF',kLDF,'qmaxN2',qmaxN2,'qmaxO2',qmaxO2, ...
                  'bN2',bN2,'bO2',bO2,'yN2_in',yN2_in, ...
                  'dp',dp,'mu',mu,'Cp_g_m',Cp_g_m,'Cp_s',Cp_s, ...
                  'lambda_eff',lambda_eff,'dH_N2',dH_N2,'dH_O2',dH_O2, ...
                  'P_in',P_in,'u_in',u_in,'T_in',T_in);

    
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.2);
    [~,Y] = ode15s(@(t,y) rhs_yN2_xT(t,y,pars), [0 tend], y0, opts);

    
    yN2_end = Y(end,1:Nx).';
    T_end   = Y(end,3*Nx+1:4*Nx).';

    prof_end = compute_pressure_velocity_x(pars, yN2_end, T_end);
    Ctot_L = prof_end.Ctot(end);          
    yO2_out = 1 - yN2_end(end);
    C_O2_out = yO2_out * Ctot_L;         
    C_N2_out = (1 - yO2_out) * Ctot_L;   
    purityC  = 100 * C_O2_out / (C_O2_out + C_N2_out);

    fprintf('Outlet O2 concentration (x=L): %.2f mol/m^3\n', C_O2_out);
    fprintf('Outlet N2 concentration (x=L): %.2f mol/m^3\n', C_N2_out);
    fprintf('Concentration-based O2 purity (x=L): %.1f%%\n', purityC);
end


function dy = rhs_yN2_xT(~,y,par)
    Nx = par.Nx; dx = par.dx; eps = par.eps; rho_p = par.rho_p;
    kLDF = par.kLDF; R = par.R;
    qmaxN2 = par.qmaxN2; qmaxO2 = par.qmaxO2; bN2 = par.bN2; bO2 = par.bO2;
    Dax = par.Dax; Cp_g_m = par.Cp_g_m; Cp_s = par.Cp_s; lambda_eff = par.lambda_eff;
    dH_N2 = par.dH_N2; dH_O2 = par.dH_O2;  yN2_in = par.yN2_in;

   
    yN2 = y(1:Nx);
    qN2 = y(Nx+1:2*Nx);
    qO2 = y(2*Nx+1:3*Nx);
    T   = y(3*Nx+1:4*Nx);

    
    prof = compute_pressure_velocity_x(par, yN2, T);
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

    
 
    Ce = eps.*Ctot.*Cp_g_m + (1-eps).*rho_p.*Cp_s;

    d2T = zeros(Nx,1);
    d2T(2:end-1) = (T(3:end) - 2*T(2:end-1) + T(1:end-2))/dx^2;
    d2T(1)   = (T(2) - T(1))/dx^2;         
    d2T(end) = (T(end-1) - T(end))/dx^2;   

    dT_dx = zeros(Nx,1);
    dT_dx(1) = (T(1) - T(1))/dx;           
    dT_dx(2:end) = (T(2:end) - T(1:end-1))/dx;

   
    Qads = (1-eps).*rho_p .* ( dH_N2.*dqN2dt + dH_O2.*dqO2dt ); 

    
    dTdt = ( lambda_eff.*d2T + Qads - (U.*Ctot*Cp_g_m).*dT_dx ) ./ Ce;

    
    dy = [dyN2dt; dqN2dt; dqO2dt; dTdt];
end


function prof = compute_pressure_velocity_x(par, yN2, T)
    Nx = par.Nx; dx = par.dx; R = par.R; eps = par.eps;
    dp = par.dp; mu = par.mu; P_in = par.P_in; u_in = par.u_in;

    
    yN2 = max(1e-9, min(1-1e-9, yN2));  yO2 = 1 - yN2;
    M_N2 = 0.0280134; M_O2 = 0.031998;
    Mmix = yN2*M_N2 + yO2*M_O2;

   
    Ctot_in = P_in/(R*T(1));
    Gm = u_in * Ctot_in;  

    P   = zeros(Nx,1);
    U   = zeros(Nx,1);
    Ctt = zeros(Nx,1);
    P(1) = P_in;

    for i = 1:Nx-1
        Ctt(i) = P(i)/(R*T(i));
        U(i)   = Gm / Ctt(i);
        rho_g  = Ctt(i) * Mmix(i);

        visc  = 150*((1-eps)^2/eps^3) * mu * U(i) / (dp^2);
        inert = 1.75*((1-eps)/eps^3) * rho_g * U(i)^2 / dp;
        dPdx  = -(visc + inert);

        P(i+1) = max(1e3, P(i) + dPdx*dx); 
    end

    
    Ctt(end) = P(end)/(R*T(end));
    U(end)   = Gm / Ctt(end);

    prof = struct('P',P,'U',U,'Ctot',Ctt);
end
