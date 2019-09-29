%% Reduced order Model for the Burgers Equation: 1D Nonlinear-Convection Diffusion %%

%close all
clear 

nu = 0.01;       % diffusion coefficient
u_init = 2;   % initial top hat level
u_l = 1;      % left boundary condition
u_r = 1;      % right boundary condition

L = 1;         % length of domain
Nx = 101;        % number of grid points
dx = L/(Nx-1);  % mesh size
tmax = 1;      % end time
Nt = 2000;       % number of timesteps
dt = tmax/(Nt);      % timestep size
u0 = [u_l*ones(1,5) u_init*ones(1,15) u_r*ones(1,81)];
%u0 = 1.0*exp(-((( dx*[0 meshgrid(1:Nx-1,1)]-0.2)/0.1).^8)) + u_l;
%u0 = [u_l*ones(1,0.5*Nx) u_r*ones(1,0.5*Nx)];

X = zeros(Nx,Nt);   %initialize Solution matrix
X(:,1) = u0;

% ===== Boundary Conditions ===== 
eastBC = 1;             % following is for periodic BC
westBC = 1;
% ===== Boundary Conditions =====

%Create the differential matrices
Dx = nu/(dx^2)*(diag(eastBC,1-Nx) + diag(westBC,Nx-1) + diag(-2*ones(Nx,1)) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1));    % x second derivative Matrix 
Dx = sparse(Dx);
Cx = -1/(2*dx)*(diag(eastBC,1-Nx) + -1*diag(westBC,Nx-1) + diag(ones(Nx-1,1),1) + diag(-1*ones(Nx-1,1),-1));            % x first derivative Matrix
Cx = sparse(Cx); 

%{
% One sided stencils for first and second derivatives
Dx(1,:) = nu/(dx^2)*[[2 -5 4 -1] zeros(1,Nx-4)];
Dx(Nx,:)= nu/(dx^2)*[zeros(1,Nx-4) [-1 4 -5 2]];
Cx(1,:) = -1/(2*dx)* [-3 4 -1 zeros(1,Nx-3)];
Cx(Nx,:) = -1/(2*dx)* [zeros(1,Nx-3) 1 -4 3];
%}

%%{
%Specification of Outflow BC:
Cx(1,:) = zeros(1,Nx);
Dx(1,:) = zeros(1,Nx);
Cx(Nx,:) = zeros(1,Nx);
Dx(Nx,:) = zeros(1,Nx);
Cx(Nx,Nx-1) = 1/dx;   
Cx(Nx,Nx) = -1/dx; 
%Cx(Nx,:) = -1/(2*dx)* [zeros(1,Nx-3) 1 -4 3];
%}

%Solve for Nt time steps
timespan = 0:dt:tmax;             % integration timespan in [seconds]
tic
[~,y] = ode45(@(t,y) odefun3FOM(y,Cx,Dx), timespan, u0);
X = y'; clear y    % solution matrix for the temporal modes a(t) at different times
toc
disp('FOM computed')

    xcoord = 1:Nx;
    figure
    plot( xcoord, X(:,2),xcoord, X(:,0.2*Nt+1), xcoord, X(:,0.4*Nt+1),xcoord, X(:,0.6*Nt+1),xcoord, X(:,0.8*Nt+1),xcoord, X(:,1*Nt+1))
    ylim([u_l*0.9,u_init*1.1])
    title('FOM Solution')
    legend('one after IC','0.4 sec','0.8 sec','1.2 sec','1.6 sec','2 sec')
    
%save('FOM_data','X','Dx','Cx','Nx','Nt','dx','dt','timespan')
    
%% POD Part %%
u0 = zeros(Nx,1);
%u0(1) = u_l;
u0 = X(:,1);
%u0 = mean(X,2);
tr = 1.0;                       % truncate the original Solution Matrix --> snapshot time factor
Xtr = X(1:Nx,1:tr*Nt+1)-repmat(u0,1,tr*Nt+1);        % Snapshot matrix up to time tr*Nt
[U,S,V] = svd(Xtr);  
[N,K] = size(Xtr);
k = 5;                         % the required k varies significantly with the amount of K snapshots! How do we estimate k with given K?

Uk = U(1:N,1:k);
Sk = S(1:k,1:k);
Vk = V(1:K,1:k);

Xlra = Uk*Sk*Vk'+repmat(u0,1,tr*Nt+1);               % low-rank-approximation of Xtr
Vlra = Sk*Vk';                     % temporal modes over snapshot time (pseudo time)
frob = norm(Xtr+repmat(u0,1,tr*Nt+1)-Xlra,'fro')

aPOD = Uk'*(X-repmat(u0,1,tr*Nt+1));      %predicted Vlra modes (exceeding snapshot time)

figure
plot(1:N, Xlra(:,0.1*Nt),'+',1:N, Xlra(:,0.4*Nt),'+',1:N, Xlra(:,0.8*Nt),'+',1:N, X(:,0.1*Nt),1:N, X(:,0.4*Nt),1:N, X(:,0.8*Nt))
title('Projection Error: Comparison Xlra to Xtr')
legend('0.1 tmax', '0.4 tmax','0.8 tmax')
%{
figure
plot(xcoord,Uk(:,1),xcoord,Uk(:,2),xcoord,Uk(:,3),xcoord,Uk(:,5),xcoord,Uk(:,7),xcoord, Uk(:,10))
title('POD Modes over x-coordinate')
legend('1st','2nd','3rd','5th','7th','10th')

tcoord = 1:K;
figure
plot(tcoord,Vlra(2,:),tcoord,Vlra(3,:),tcoord,Vlra(4,:),tcoord,Vlra(5,:),tcoord,Vlra(10,:))
title('Vlra Modes up to cut off time')
legend('2nd','3rd','4th','5th','10th')
%}

save('FOM_data','X','Dx','Cx','Nx','Nt','dx','dt','timespan','Uk','Sk','Vk')

%% Precomputation & Galerkin Projection %%
A0star = Dx + diag(u0)*Cx;
Dxrom = Uk'*Dx*Uk;            % GALERKIN PROJECTION - obtain the kxk matrix to solve for arom(1) to arom(k)
CxUk = Cx*Uk;

% Precompute QX advective matrix for nonlinear terms
qX = zeros(Nx,k^2);
for ii = 1:k
    for jj = 1:k
       qX(:,(ii-1)*k+jj) = Uk(:,ii).* CxUk(:,jj);  %non conservative 
       %qX(:,(ii-1)*k+jj) = Uk(:,ii).* Uk(:,jj);     %conservative 
    end
end
clear ii jj
QX = Uk' *qX;          %non conservative 
%QX = 0.5*Uk' *Cx *qX;   %conservative 
Qrom = QX; clear qX           % Completely precomputed quadratic matrix for nonlinear terms

% Precompute other terms due to constant subtraction
y0 = u0;
Cxy0 = Cx*y0;
D0 = Uk'*Dx*y0;

% LX2
lX2 = zeros(Nx,k);
for ii = 1:k
    lX2(:,ii) = Uk(:,ii).* Cxy0(:,1);  %non conservative
    %lX2(:,ii) = Uk(:,ii).* y0(:,1);     %conservative
end
clear ii
LX2 = Uk'* lX2;            %non conservative 
%LX2 = 0.5 *Uk' *Cx *lX2 ;  %conservative --> split in half for formal reasons
L2 = LX2; clear lX2

% LX3
lX3 = zeros(Nx,k);
for ii = 1:k
    lX3(:,ii) = y0(:,1).* CxUk(:,ii);
end
clear ii
LX3 = Uk'* lX3;            %non conservative
%LX3 = LX2;                  %conservative --> other half of LX2
L3 = LX3; 
Lrom = L2 + L3 + Dxrom; clear lX3 L2 L3

% LX4
LX4 = Uk'* (y0.* Cxy0);        %non conservative
%LX4 = 0.5 *Uk' *Cx *(y0.*y0);   %conservative
Crom = LX4 + D0; clear LX4   

%% Standard Galerkin ROM %%
dt_int = dt;                          % choose time step size for the integration to get arom
Na_t = tmax/dt_int;
timespan = [0:dt_int:tmax];             % integration timespan in [seconds]
a0_ROM = Uk'*Xtr(:,1);
% [t,y] = ode45(@(t,y) odefun3(y,A0star,Uk,Cx,u0,Qrom,Lrom,Crom), timespan, a0_ROM);
% arom = y';    % solution matrix for the temporal modes a(t) at different times
arom = zeros(k,length(timespan));
arom(:,1)= a0_ROM;
for tt = 1:length(timespan)-1
    arom(:,tt+1) = arom(:,tt) + dt*(Qrom * kron(arom(:,tt),arom(:,tt)) + Lrom*arom(:,tt) + Crom);
end
u_rom = Uk*arom + repmat(y0,1,Na_t + 1);

save('temporal_modes','aPOD','arom');

%% Results Plotting Part %%
time_predict1 = 0.25;       % solution prediction time 1
time_predict2 = 0.5;       % solution prediction time 2
time_predict3 = 1.0;       % solution prediction time 3
tt1 =round(time_predict1*Na_t/tmax)+1;
tt2 =round(time_predict2*Na_t/tmax)+1;
tt3 =round(time_predict3*Na_t/tmax)+1;

figure
plot(xcoord,u_rom(:,tt1),'-o',xcoord,X(:,tt1),xcoord,u_rom(:,tt2),'-o',xcoord, X(:,tt2),xcoord,u_rom(:,tt3),'-o',xcoord,X(:,tt3));  % compare ROM prediction with actual full solution at different (future) times
title('Predicted Solution Vs Full Solution, variable space, fixed time')
legend('ROM Solution','FOM Solution')

tt = 1:length(X);
figure
plot(tt,aPOD(2,:),'-.',tt,aPOD(3,:),'-.',tt,aPOD(4,:),'-.',tt,aPOD(5,:),'-.',...
    tt,arom(2,:),tt,arom(3,:),tt,arom(4,:),tt,arom(5,:))
title('Predicted Vlra modes vs a(t) (exceeding snapshot time)')
legend('2nd aPOD','3rd aPOD','4th aPOD','5th aPOD','2nd a(t)','3rd a(t)','4th a(t)','5th a(t)')

%{
figure
pseudo_dt = round(1:tr*Na_t);
plot(pseudo_dt,arom(2,pseudo_dt),pseudo_dt,arom(3,pseudo_dt),pseudo_dt,arom(4,pseudo_dt),pseudo_dt,arom(5,pseudo_dt),pseudo_dt,arom(10,1:tr*Na_t))
title('a(t) temporal modes over 3 seconds')
legend('2nd','3rd','4th','5th','10th')
%}

%% Eddy Viscosity %% 
Dt = 1/(2*dt)*(diag(ones(Nt-1,1),1) + diag(-1*ones(Nt-1,1),-1));            % x first derivative Matrix
Dt = sparse(Dt); %Dt(1,1)= -1/(2*dt); Dt(Nt,Nt)= 1/(2*dt); 
Dt(1,:) = 1/(2*dt)* [-3 4 -1 zeros(1,Nt-3)];
Dt(Nt,:) = 1/(2*dt)* [zeros(1,Nt-3) 1 -4 3];
adotPOD = (Dt*aPOD(:,1:end-1)')';
adotrom = (Dt*arom(:,1:end-1)')';

%adotPOD = Uk'*(Dx*X + X.*(Cx*X)); adotPOD(:,end-1) = [];

%adotPOD = zeros(size(aPOD(:,1:end-1)));
%adotrom = zeros(size(aPOD(:,1:end-1)));
for tt = 1:length(adotrom)
%     adotPOD(:,tt) = (aPOD(:,tt+1) - aPOD(:,tt))/dt;
%     adotrom(:,tt) = (arom(:,tt+1) - arom(:,tt))/dt;
%     adotPOD(:,tt) = Qrom * kron(aPOD(:,tt),aPOD(:,tt)) + Lrom*aPOD(:,tt) + Crom;
%     adotrom(:,tt) = Qrom * kron(arom(:,tt),arom(:,tt)) + Lrom*arom(:,tt) + Crom;
end
tt = 1:length(adotrom);
figure
plot(tt,adotPOD(1,:),tt,adotPOD(2,:),tt,adotPOD(3,:),tt,adotPOD(4,:),...
    tt,adotrom(1,:),'-.',tt,adotrom(2,:),'-.',tt,adotrom(3,:),'-.',tt,adotrom(4,:),'-.')
title('d/dt(a_{POD}) vs d/dt(a_{ROM}) discrepancy')

%figure
% nu_t = zeros(size(arom));
% for tt = 1:length(adotrom)
%     nu_t(:,tt) = (adotPOD(:,tt) - adotrom(:,tt)).*((Dxrom*arom(:,tt) + Crom).^(-1)); %  (Dxrom/nu *arom(:,tt)).*nu_t(:,tt))
      %plot(1:tt,nu_t(1,1:tt),1:tt,nu_t(2,1:tt),1:tt,nu_t(3,1:tt),1:tt,nu_t(4,1:tt),1:tt,nu_t(5,1:tt))
% end

nu_t = (adotPOD - adotrom).*((abs(Dxrom*arom(:,1:end-1))/nu + abs(D0/nu)).^(-1));
max(max(nu_t))
tmp1 = (abs(Dxrom*arom(:,1:end-1))/nu + abs(D0)/nu).^(1);
tmp1 = (Dxrom*arom(:,1:end-1)/nu + D0/nu).^(1);
figure; grid on; 
plot(tt,tmp1(1,:),tt,tmp1(2,:),tt,tmp1(3,:),tt,tmp1(4,:),tt,tmp1(5,:))
title('Linear + Const term in denominator')

figure; 
plot(tt,nu_t(1,:),tt,nu_t(2,:),tt,nu_t(3,:),tt,nu_t(4,:),tt,nu_t(5,:))
title('Modal "artificial" viscosity')

%% ROM with eddy viscosity "closure"%%
arom_t = zeros(size(arom));
arom_t(:,1)= a0_ROM;
for tt = 1:length(adotrom)
    %arom_t(:,tt+1) = arom(:,tt) + dt*(Qrom * kron(aPOD(:,tt),aPOD(:,tt)) + Lrom*aPOD(:,tt) + Crom);
    arom_t(:,tt+1) = arom_t(:,tt) + dt*(Qrom * kron(arom(:,tt),arom(:,tt)) + Lrom*arom(:,tt) + Crom... + nu_t(:,tt));
                    + ( abs(Dxrom*arom(:,tt)/nu) + abs(D0/nu) ).*nu_t(:,tt)); %  + Crom/nu.*nu_t(:,tt)
end 
max(max(arom_t))

tt = 1:length(X);
figure
plot(tt,aPOD(2,:),'-.',tt,aPOD(3,:),'-.',tt,aPOD(4,:),'-.',tt,aPOD(5,:),'-.',...
    tt,arom_t(2,:),tt,arom_t(3,:),tt,arom_t(4,:),tt,arom_t(5,:))
title('Artificial visc closure ROM')
legend('2nd aPOD','3rd aPOD','4th aPOD','5th aPOD','2nd a(t)','3rd a(t)','4th a(t)','5th a(t)')