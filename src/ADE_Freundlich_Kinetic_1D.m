% ADE_Freundlich_Kinetic.m
% date: 22 April 2025
% author: Zijie and Marc

% Description: This function is used for gridsearch to fit the experiment
% data, in order to find the best f and alpha value 

% INPUT:
%   Q      - flow rate [mL/h]
%   Kf     - Freundlich coefficient [(g/kg)/(g/m^3)^(nf)]
%   nf     - Freundlich exponent [-]
%   alpha  - mass transfer rate constant [1/s or similar units]
%   f      - fraction of equilibrium sorption sites [-]
%   rho    - bulk density [kg/m^3]
%
% OUTPUT:
%   time   - time vector [s]
%   c_out  - concentration at outlet vs. time [g/m^3]


function [time, c_out] = ADE_Freundlich_Kinetic_1D(c0, Q, rho, Kf, nf, f, alpha,t_inject,tmax)
% Conversions
hr2s = 60^2;
dy2s  = 24*hr2s;
wk2s  = 7*dy2s;

% Physical Paramters
Length = 0.07; % [m] Column length
theta = 0.33;  % [-] Porosity
A = 1.766e-4;        % [m^2]
qb = (Q/1e6 / hr2s)*dy2s / A;          % Darcy flux
v = qb / theta;      % interstitial velocity
Disp = 0.006 * Q/12; % [m^2/s] dispersion
cini = 1e-6; 

smax = (1-f)*Kf*c0^nf;
sini = (1-f)*Kf*cini^nf;

% Numerical parameters
Nt = 400;
gamma = 0;  % 1 = FE (explicit), 0 = BE (implicit), 0.5 = C-N (implicit)
dt = tmax/Nt;
epsilon = 1e-8;
tol = 1e-8;
m_max = 15;

% Peclet number
Pe = v*Length/(theta*Disp);
%fprintf('Peclet number: Pe = %3.2f.\n',Pe)

% Build Grid
Grid.xmin = 0; Grid.xmax = Length; Grid.Nx = 1e3;
Grid = build_grid(Grid);
dt_cfl = Grid.dx/v;
if gamma == 1 && dt > dt_cfl
    error('dt exceeds maximum timestep!')
else
    fprintf('Timestep %3.2f dt_cfl.\n\n',dt/dt_cfl)
end

% Build operatore
[Div,Grad,~,I,M] = build_ops(Grid);
fs = zeros(Grid.Nx,1);
q = v*ones(Grid.Nf,1);
A = flux_upwind(q,Grid);
Theta = theta*I;
Kd = Disp*speye(Grid.Nf);
L = theta*Div*(A - Kd*Grad);

% Boundary conditions
BC.dof_dir = Grid.dof_xmin;
BC.dof_f_dir = Grid.dof_f_xmin;
BC.g = 0;  % This is for the N-R update


[B,N,fn,BC] = build_bnd(BC,Grid,I);
% Modify to accounf for 's' unknowns
B = [B 0*B];
N = [N,spalloc(Grid.Nx,Grid.Nx,0);...
    spalloc(Grid.Nx,Grid.Nx-1,0),I];

%% Solve time evolution
c = cini+zeros(Grid.N,1); cold = c; c(1)=c0;
s = sini*ones(Grid.N,1); sold = s; %s(1) = smax;
c = cini + cini/10*sin(2*pi*Grid.xc/Length); c(1)=c0; 


u = [c;s]; uold = [cold;sold];
res = residual_ade_freundlich_kinetic(u,uold,L,Theta,fs,dt,gamma,rho,Kf,nf,alpha,f,Grid);
Jac_ana = jacobian__ade_freundlich_kinetic_ana(u, uold, L, Theta, fs, dt, gamma, rho, Kf, nf, alpha, f, Grid);
Jac_num = jacobian__ade_freundlich_kinetic_num(u, uold, L, Theta, fs, dt, gamma, rho, Kf, nf, alpha, f, Grid, epsilon);


c_x07 = zeros(Nt+1, 1);
x_target = 0.07;
[~, index_x] = min(abs(Grid.xc - x_target));


c_out = cini+zeros(Nt+1,1);  % effluent concentration

u = [c;s]; 
for n = 1:Nt
     t = (n-1)*dt; % current time
    if t <= t_inject
       u(1) = c0;
    else
       u(1) = 0;
    end
 
    % Newton-Raphson Iteration
    nres = 1; ndu = 1; m = 0;
    uold = u;
    res = residual_ade_freundlich_kinetic(u,uold,L,Theta,fs,dt,gamma,rho,Kf,nf,alpha,f,Grid);
    while (nres > tol || ndu > tol) && m < m_max

        Jac_ana = jacobian__ade_freundlich_kinetic_ana(u, uold, L, Theta, fs, dt, gamma, rho, Kf, nf, alpha, f, Grid);        

        du = solve_lbvp(Jac_ana,-res,B,BC.g,N);
        u = u + du; 
        res = residual_ade_freundlich_kinetic(u,uold,L,Theta,fs,dt,gamma,rho,Kf,nf,alpha,f,Grid);
        nres = norm(N'*res); ndu = norm(N'*du);
        m = m+1;
%        fprintf('it = %d: time = %3.2f: nres = %3.2e  ndc = %3.2e\n',m,n*dt,nres,ndu)
        if m == 1; ndtheta = 0; end % to allow exit on first iteration
    end
    if m<m_max;
%        fprintf('Newton converged after %d iterations.\n\n',m)
    else
        fprintf('Newton did not converge\n\n')
    end

    if mod(n,10) == 0
        c = u(1:Grid.Nx); s = u(Grid.Nx+1:2*Grid.Nx);

    end
    c_out(n+1) = u(Grid.Nx);
    c_x07(n+1) = u(index_x);
end


 time = 0:dt:tmax;
 time_hrs = time/hr2s;
 time_dy = time/dy2s;


%% Residual
function [res] = residual_ade_freundlich_kinetic(u,uold,L,Theta,fs,dt,gamma,rho,Kf,nf,alpha,f,Grid)
% Residual for the ADE with both Freundlich and kinetic sorption
% 0) Extrace variables
c = u(1:Grid.Nx);               cold = uold(1:Grid.Nx); 
s = u(Grid.Nx+1:2*Grid.Nx);     sold = uold(Grid.Nx+1:2*Grid.Nx);

% 1) Conservation of solute c
gam_c = gamma*cold+(1-gamma)*c;
gam_s = gamma*sold+(1-gamma)*s;
Theta_inv = inv(Theta);

res_c = Theta*c + f*rho*Kf*c.^nf - (Theta*cold + f*rho*Kf*cold.^nf) ...
      + dt*L*gam_c - dt*fs + dt*alpha*rho*( (1-f)*Kf*gam_c.^nf-gam_s );


res_c(1) = 0; % hack! (BC)

% 2) conservation of sorbed s
res_s = s - sold - dt*alpha*((1-f)*Kf*gam_c.^nf-gam_s);

% 3) Assemble full residual
res = [res_c;res_s];
end


%% Jacobian Analytical
function [Jac_ana] = jacobian__ade_freundlich_kinetic_ana(u,uold,L,Theta,fs,dt,gamma,rho,Kf,nf,alpha,f,Grid)

% Extract variables
c = u(1:Grid.Nx);               
cold = uold(1:Grid.Nx); 
s = u(Grid.Nx+1:2*Grid.Nx);     
sold = uold(Grid.Nx+1:2*Grid.Nx);

% Compute intermediate terms
gam_c = gamma*cold + (1-gamma)*c;
gam_s = gamma*sold + (1-gamma)*s;
Theta_inv = inv(Theta);

% Compute derivatives
dc_eqbm = nf * f * rho * Kf * c.^(nf - 1);
dc_kin_c = nf * (1-f) * Kf * (1-gamma) * gam_c.^(nf-1);

Jcc = Theta + spdiags(dc_eqbm, 0, Grid.Nx, Grid.Nx) + dt * (1 - gamma) * L ...
       + dt * alpha * rho *spdiags(dc_kin_c, 0, Grid.Nx, Grid.Nx);
Jcs = -dt * alpha * rho * (1-gamma)*speye(Grid.Nx);  % Partial derivative of c equation w.r.t. s
Jsc = -dt * alpha * spdiags(dc_kin_c, 0, Grid.Nx, Grid.Nx);  % Partial derivative of s equation w.r.t. c
Jss = speye(Grid.Nx) + dt * alpha * (1-gamma) * speye(Grid.Nx); % Partial derivative of s equation w.r.t. s

% Assemble full Jacobian
Jac_ana = [Jcc, Jcs;
            Jsc, Jss];

end
%end

% % %% Jacobian Numerical
function [Jac_num] = jacobian__ade_freundlich_kinetic_num(u,uold,L,Theta,fs,dt,gamma,rho,Kf,nf,alpha,f,Grid,epsilon)
u_perturb=u;

% Pre-allocate storage (important to speed up the for loop)
Jac_num = spalloc(2*Grid.N,2*Grid.N,3*(2*Grid.N));
res = zeros(2*Grid.N,1);
res_perturb = res;
% 
% %% Computing the Jacobian by finite difference
% % Loop over each unknown and
% % 1) Perturb it
% % 2) Compute the change in residual vector
% % 3) Store the change from unperturbed vector as column of Jacobnian
% % 4) Reset the perturbed value
% 
for i=1:2*Grid.N
    u_perturb(i)=u(i)+epsilon;
    res         = residual_ade_freundlich_kinetic(u        ,uold,L,Theta,fs,dt,gamma,rho,Kf,nf,alpha,f,Grid);
    res_perturb = residual_ade_freundlich_kinetic(u_perturb,uold,L,Theta,fs,dt,gamma,rho,Kf,nf,alpha,f,Grid);
    Jac_num(:,i)=(res_perturb-res)/epsilon;
    u_perturb(i)=u(i);
end
end
end