%% ROM Calibration / Closure Modeling routine for 2D Navier Stokes equation/ Burgers equation%%

clear
%close all

load('uniformsnapshots_128_Re30000_100sec.mat')
load('FOM_data2D_128.mat')
X = [XU_uni;XV_uni];
Xmean = mean(X,2);
M_size = 2*Nx*Ny;

%% POD Part %%
Nt = Nt-1;
X(:,1) = [];
tr = 1.0;                       % Truncation time factor - trucates original Solution Matrix X
y0 = X(:,1);          % Add or remove percent sign in order to suppress or enable subtraction of IC 
trvec = 10:10:tr*Nt;    % Vector of indeces to be selected from the snapshot matrix
Xtr = X(:,trvec)-repmat(y0,1,length(trvec));        % Snapshot matrix up to time tr*Nt+1, note: first column is initial condition
tic
[U,S,V] = svd(Xtr,'econ'); 
disp('SVD computed')
toc
[N,K] = size(Xtr);
k = 50;                         % Number of chosen POD modes

Uk = U(1:N,1:k);
Sk = S(1:k,1:k);
Vk = V(1:K,1:k);

Xlra = Uk*Sk*Vk'+repmat(y0,1,length(trvec));               % Low-rank-approximation of Xtr
Vlra = Sk*Vk';                  % temporal modes over snapshot time (pseudo time)
frob = norm(Xtr+repmat(y0,1,length(trvec))-Xlra,'fro')
POD_energy = (trace(S(1:min(N,K),1:min(N,K)))-trace(Sk))/trace(S(1:min(N,K),1:min(N,K)))

aPOD = Uk'*(X-repmat(y0,1,size(X,2)));                  %predicted Vlra modes (exceeding snapshot time)
Xlrau = Xlra(1:Nx*Ny,:);
Xlrav = Xlra((Nx*Ny+1):(2*Nx*Ny),:);

TKE_POD = sum(aPOD.^2);

%{
figure
plot(1:size(TKE_POD,2),TKE_POD)
title('TKE_{POD}')

figure
plot(1:K,Vlra(1,:),1:K,Vlra(2,:),1:K,Vlra(3,:),1:K,Vlra(4,:),1:K,Vlra(5,:))
title('Temporal Modes up to cut off time')
legend('1st','2nd','3rd','4th','5th')

figure
for tt = 1:0.1*(K-1):K
    %%{
    p = quiver(reshape(Xlrau(:,tt),[Ny,Nx]),reshape(Xlrav(:,tt),[Ny,Nx]),3); colorbar; shading interp; colormap('jet'); axis equal; set(gca,'Ydir','reverse');
    drawnow
    %{
    basename1 = 'ProjErrVec';
    filename1 = [basename1,num2str(tt),'.jpg'];
    saveas(p,filename1);
    %}
end
clear tt
%}

%% ROM/Galerkin Projection Part %%

Dstar = Uk'*D*Uk;            % GALERKIN PROJECTION - obtain the kxk matrix to solve for a_t(1) to a_t(k)
Uk1 = Uk(1:Nx*Ny,:);
Uk2 = Uk((Nx*Ny+1):(2*Nx*Ny),:);
Uk11 = [Uk1;Uk1];
Uk22 = [Uk2;Uk2];
CXUk = CX*Uk;
CYUk = CY*Uk;

% Precompute QX and QY advective matrices for nonlinear terms
qX = zeros(M_size,k^2);
qY = zeros(M_size,k^2);
for ii = 1:k
    for jj = 1:k
       qX(:,(ii-1)*k+jj) = Uk11(:,ii).* CXUk(:,jj);
       qY(:,(ii-1)*k+jj) = Uk22(:,ii).* CYUk(:,jj);
    end
end
clear ii jj
QX = Uk' *qX;
QY = Uk' *qY;
Qrom = QX + QY;            % Completely precomputed quadratic matrix for nonlinear terms

Q3D = zeros(k,k,k);
for m = 1:k
    Q3D(:,:,m) = Qrom(:, ((m-1)*k+1) : m*k);
end
clear m

% Precompute other terms due to constant subtraction
D0 = Uk'*D*y0;
u0 = y0(1:Nx*Ny,1); 
v0 = y0((Nx*Ny+1):(2*Nx*Ny),:);
uu0 = [u0; u0];
vv0 = [v0; v0];
CXuv = CX*y0;
CYuv = CY*y0;

% LX2 & LY2
lX2 = zeros(M_size,k);
lY2 = zeros(M_size,k);
for ii = 1:k
    lX2(:,ii) = Uk11(:,ii).* CXuv(:,1);
    lY2(:,ii) = Uk22(:,ii).* CYuv(:,1);
end
clear ii
LX2 = Uk'* lX2;
LY2 = Uk'* lY2;
L2 = LX2 + LY2;

% LX3 & LY3
lX3 = zeros(M_size,k);
lY3 = zeros(M_size,k);
for ii = 1:k
    lX3(:,ii) = uu0(:,1).* CXUk(:,ii);
    lY3(:,ii) = vv0(:,1).* CYUk(:,ii);
end
clear ii
LX3 = Uk'* lX3;
LY3 = Uk'* lY3;
L3 = LX3 + LY3;

% LX4 & LY4 
LX4 = Uk'* (uu0.* CXuv);
LY4 = Uk'* (vv0.* CYuv);
Crom = LX4 + LY4 + D0;

Lrom = L2 + L3 + Dstar;          % Completely precomputed linear matrix for linear terms 

t_max = 100.0;                              % total simulation time
dt_int = 0.01;                          % integration time size (not actual time discretization step size)
Na_t = t_max/dt_int;                      % number of temporal steps/ snapshots! 
timespan = dt:dt_int:t_max;             % integration timespan in [seconds]     

%% Optimization Part %%

%{
x0 = zeros(k+k^2+k^3,1);
lb = -1*abs(max([Crom(:) ; Lrom(:) ; Qrom(:)])*ones(size(x0)));
ub = 1*abs(max([Crom(:) ; Lrom(:) ; Qrom(:)])*ones(size(x0)));
%}
%{
x0 = zeros(k+k^2,1);
lb = -1*abs(max([Crom(:) ; Lrom(:)])*ones(size(x0)));         % Lower bound
ub = 1*abs(max([Crom(:) ; Lrom(:)])*ones(size(x0)));          % Upper bound
%}
%%{
%x0 = zeros(6*k,1);
%x0 = zeros(k+2*k+(k+k^2),1);
%x0 = zeros((2*k+2*k^2),1);
%x0 = [0.01*max(Crom)*rand(k,1) ; 0.01*max(max(Lrom))*rand(2*k,1)];
%x0 = [0.001*max(Crom)*rand(k,1) ; 0.001*max(max(Lrom))*rand(2*k,1) ; 0.00001*max(max(Qrom))*rand(3*k,1)];
x0 = [0.01*max(Crom)*rand(k,1) ; 0.0001*max(max(Lrom))*rand(2*k,1) ; 0.00001*max(max(Qrom))*rand(3*k,1)];
%x0 = x + 0.000001*max(max(Qrom))* rand(size(x));
%x0 = [0.0001*max(Crom)*rand(k,1) ; 0.0001*max(max(Lrom))*rand(k^2,1) ; 0.0001*max(max(Qrom))*rand(k+k^2,1)];
lb = -1*abs(max([Crom(:) ; Lrom(:)])*ones(size(x0)));         % Lower bound
ub = 1*abs(max([Crom(:) ; Lrom(:)])*ones(size(x0)));          % Upper bound
%}
a0_ROM = Uk'*Xtr(:,1);
m = 30;

%%{ 
options = optimoptions(@lsqnonlin,'Display','iter');
options.Algorithm = 'levenberg-marquardt';
options.MaxFunEvals = 50*size(x0,1);
options.UseParallel = true;
options.InitDamping = 1.5;
f = @(x)objfun(x,Qrom,Lrom,Crom,k,timespan,a0_ROM,aPOD,TKE_POD,Uk,y0,X,Nx,Q3D,m);
[x,resnorm,residual,exitflag,output] = lsqnonlin(f,x0,[],[],options);
output.firstorderopt

%standard full rank Q,L,C
%{
Q_final = reshape(x((k+k^2+1):(k+k^2+k^3),1),[k,k^2]);
%Q_final = zeros(k,k^2);                         %Q optimization disabled
L_final = reshape(x((k+1):(k+k^2),1),[k,k]);
C_final = reshape(x(1:k,1),[k,1]);
%}

%low rank Q & L
%%{
%Q_final = zeros(k,k^2);
%%{
Q_final = 0;
L_final = x((k+1):(2*k),1)*x((2*k+1):(3*k),1)';
C_final = x(1:k,1);
%Q3D_final = zeros(k,k,k);
Q{1} = x((3*k+1):(4*k),1);
Q{2} = x((4*k+1):(5*k),1);
Q{3} = x((5*k+1):(6*k),1);
Q3D_final = cpdgen(Q);
%}
%}

%low rank Q only
%{
Q_final = x((k+k^2+1):(k+k^2+k))*x((2*k+k^2+1):(2*k+2*k^2))';
L_final = reshape(x((k+1):(k+k^2)),[k,k]);
C_final = x(1:k);
%}

%%{
[t,a_dot] = ode45(@(t,a) odefun3calibrate(a,Qrom,Lrom,Crom,Q_final,L_final,C_final,Q3D,Q3D_final,k), timespan, a0_ROM);
aROM_opti = a_dot';    % solution matrix for the temporal modes a(t) at different times
%}
Na_t = Nt;
tmax = dt*Nt;

f_rom = Uk*aROM_opti + repmat(y0,1,size(t,1));
u_rom = f_rom(1:Nx*Ny,:);
v_rom = f_rom((Nx*Ny+1):(2*Nx*Ny),:);

save('calibrateA_Re30000_lowrankQLC_k25','Q3D_final','L_final','C_final','aROM_opti','aPOD','u_rom','v_rom','x')

%% Results Plotting Part %%

TKE_ROMopti = sum(aROM_opti.^2);

tt = 1:length(aPOD);
figure
plot(tt,TKE_POD,tt,TKE_ROMopti,'o')
title('TKE_{POD} vs TKE_{ROM}')
legend('TKE_{POD}','TKE_{ROM,optimal}')
clear tt

%load('VelocityROM_U_V.mat')
for tt = 1:round(0.05*Na_t):Na_t+1
    %%{
    p = quiver(reshape(u_rom(:,tt),[Ny,Nx]),reshape(v_rom(:,tt),[Ny,Nx]),3); colorbar; shading interp; colormap('jet'); axis equal; set(gca,'Ydir','reverse');
    drawnow
    %{
    basename1 = 'LDC_ROM_vec';
    filename1 = [basename1,num2str(tt),'.jpg'];
    saveas(p,filename1);
    %}
end
clear tt

%%{
nu = 1/30000;
Dx_int = 1/nu*Dx(2:Nx-1,2:Nx-1);
Dy_int = 1/nu*Dy(2:Ny-1,2:Ny-1);
tic
tmp = full(kron(Dx_int,eye(Ny-2)) + kron(eye(Nx-2),Dy_int));
D_psi_int = (tmp)\eye((Nx-2)*(Ny-2));
disp('inverse laplacian computed')      %needs 200 seconds for a 128x128 grid
toc
clear tmp
Psi_rom = zeros(Ny,Nx);
Psi_rom_int = Psi_rom(2:Ny-1,2:Nx-1);
%}

for tt = 1:round(0.05*Na_t):Na_t+1
    %%{
    U_rom = reshape(u_rom(:,tt),[Ny,Nx]);
    V_rom = reshape(v_rom(:,tt),[Ny,Nx]);
    W_rom = -V_rom*Cx'+ Cy*U_rom;
    W_rom_int = W_rom(2:Ny-1,2:Nx-1);
    Psi_rom_int(:) = -D_psi_int*W_rom_int(:);   % compute psi on interior domain only through in inverse laplacian
    Psi_rom(2:Ny-1,2:Nx-1) = Psi_rom_int;
    
    p = pcolor(Psi_rom); colorbar; shading interp; colormap('jet'); axis equal; set(gca,'Ydir','reverse');
    drawnow
    %{
    basename1 = 'LDC_ROM_psi';
    filename1 = [basename1,num2str(tt),'.png'];
    saveas(p,filename1);
    %}
end
clear tt

tt = 1:size(aPOD,2);
figure
plot(tt,aPOD(1,tt),'b-.',tt,aPOD(2,tt),'r-.',tt,aPOD(3,tt),'g-.',tt,aPOD(4,tt),'c-.',tt,aPOD(5,tt),'y-.',tt,aROM_opti(1,tt),'b',tt,aROM_opti(2,tt),'r',tt,aROM_opti(3,tt),'g',tt,aROM_opti(4,tt),'c',tt,aROM_opti(5,tt),'y')
title('Calibrated a_{ROM} vs a_{POD} LDC'); xlabel('time'); ylabel('a(t)'); 
legend('1st POD','2nd POD','3rd POD','4th POD','5th POD','1st ROM','2nd ROM','3rd ROM','4th ROM','5th ROM')