%% Reduced Order Model for 2D Burgers Equation with BETTER Data Management %%
clear
close all 

%% Pre Setup %%
nu = 1.0/1000;             

L_x = 1;       % Domain Length in x-Direction
L_y = 1;        % Domain Length in y-Direction
Nx = 100;      
Ny = 100;
dx = L_x/(Nx-1)
dy = L_y/(Ny-1)

%% Creation of system matrices for FOM %%
% ===== Boundary Conditions ===== 
%%{
eastBC = zeros(Ny-2,1);             % for Dirichlet BC's
westBC = zeros(Ny-2,1);
northBC = 3.0*exp(-((( dx*[0 meshgrid(1:Nx-1,1)]-0.5)/0.4).^16)); northBC = northBC(2:end-1)'; 
southBC = zeros(Nx-2,1);
%}
% ===== Boundary Conditions =====

Dx = nu/dx^2*(diag(-2*ones(Nx,1)) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1));                  % x second derivative Matrix
Dx = sparse(Dx);
Dy = nu/dy^2*(diag(-2*ones(Ny,1)) + diag(ones(Ny-1,1),1) + diag(ones(Ny-1,1),-1));  % y second derivative Matrix
Dy = sparse(Dy);

c_x = 1; c_y = 1;
Cx = c_x/(2*dx)*(diag(-1*ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1));         % x first derivative Matrix
Cx = sparse(Cx); 
Cy = c_y/(2*dy)*(diag(-1*ones(Ny-1,1),1) + diag(ones(Ny-1,1),-1));         % y first derivative Matrix
Cy = sparse(Cy);

Dx_int = Dx(2:Nx-1,2:Nx-1);
Dy_int = Dy(2:Ny-1,2:Ny-1);
BDxEW = nu/(dx^2)*[westBC zeros(Ny-2,Nx-4) eastBC];
BDyNS = nu/(dy^2)*[northBC'; zeros(Ny-4,Nx-2); southBC'];

Cx_int = Cx(2:Nx-1,2:Nx-1);
Cy_int = Cy(2:Ny-1,2:Ny-1);
BCxEW = c_x/(2*dx)*[westBC zeros(Ny-2,Nx-4) -eastBC];
BCyNS = c_y/(2*dy)*[northBC'; zeros(Ny-4,Nx-2); -southBC'];

%%{
Dx(1,:)= zeros(1,Nx);   Dx(Nx,:)= zeros(1,Nx);
Dy(1,:)= zeros(1,Ny);   Dy(Ny,:)= zeros(1,Ny);
Cx(1,:)= zeros(1,Nx);   Cx(Nx,:)= zeros(1,Nx);
Cy(1,:)= zeros(1,Ny);   Cy(Ny,:)= zeros(1,Ny);
%}

%{
Dx(1,:) = nu/(dx^2)*[[1 -2 1] zeros(1,Nx-3)];
Dx(Nx,:)= nu/(dx^2)*[zeros(1,Nx-3) [1 -2 1]];
Dy(1,:) = nu/(dy^2)*[[1 -2 1] zeros(1,Ny-3)];
Dy(Ny,:)= nu/(dy^2)*[zeros(1,Ny-3) [1 -2 1]];
Cx(1,:) = -1/(2*dx)* [-1 1 zeros(1,Nx-2)];
Cx(Nx,:) = -1/(2*dx)* [zeros(1,Nx-2) -1 1];
Cy(1,:) = -1/(2*dy)* [-1 1 zeros(1,Ny-2)];
Cy(Ny,:) = -1/(2*dy)* [zeros(1,Ny-2) -1 1];
%}

%{
Dx(1,:) = nu/(dx^2)*[[2 -5 4 -1] zeros(1,Nx-4)];  Dy(1,:) = nu/(dy^2)*[[2 -5 4 -1] zeros(1,Ny-4)];
Dx(Nx,:)= nu/(dx^2)*[zeros(1,Nx-4) [-1 4 -5 2]];  Dy(Ny,:)= nu/(dy^2)*[zeros(1,Ny-4) [-1 4 -5 2]];

Cx(1,:) = -1/(2*dx)* [-3 4 -1 zeros(1,Nx-3)]; Cy(1,:) = -1/(2*dy)* [-3 4 -1 zeros(1,Ny-3)];
Cx(Nx,:) = -1/(2*dx)* [zeros(1,Nx-3) 1 -4 3]; Cy(Ny,:) = -1/(2*dy)* [zeros(1,Ny-3) 1 -4 3];
%}

%% Solving the FOM %%
%{
u0 = zeros(Ny,Nx);                   % initial condition
v0 = zeros(Ny,Nx);
for jj = 1:Ny      
    for ii = 1:Nx
        xi = ii/Nx*L_x;
        yj = jj/Ny*L_y;
        
        %u0(jj,ii) = exp(-(((xi-2.5)/0.5)^8 + ((yj-2.5)/0.5)^8)) + 1.5;
        %v0(jj,ii) = u0(jj,ii);
        
        u0(jj,ii) = exp(-(((xi-2.5)/0.5)^8 + ((yj-2.5)/0.5)^8)) + 1.5;
        v0(jj,ii) = exp(-(((xi-2.5)/0.5)^8 + ((yj-2.5)/0.5)^8)) + 1.5;
        
    end
end
figure
pcolor(u0)            % test the initial condition
axis equal
set(gca,'Ydir','reverse')
%surf(y0)
clear ii jj
%}
t_max = 5.0;                              % total simulation time
dt = 0.001;                          % integration time size (not actual time discretization step size)
Nt = t_max/dt;                      % number of temporal steps/ snapshots! 
timespan = 0:dt:t_max;             % integration timespan in [seconds]

% Assemble BIG Block Matrices for System
CX = blkdiag(kron(Cx,eye(Ny)),kron(Cx,eye(Ny)));
CY = blkdiag(kron(eye(Nx),Cy),kron(eye(Nx),Cy));
D = blkdiag(kron(Dx,eye(Ny)) + kron(eye(Nx),Dy) , kron(Dx,eye(Ny)) + kron(eye(Nx),Dy));

%{
tic
f0 = [u0(:);v0(:)];
[t,f] = ode45(@(t,f) ode2DBurgers(f,CX,CY,D,Nx,Ny), timespan, f0);  %solve FOM 2D Burgers equation
X = f';                                                             % solution matrix ==> for snapshot matrix
%{
u = u0;
Xu = u0(:);
v = v0;
Xv = v0(:);
for nn = 1:Nt
    u = u + dt*(u.*(u*Cx') + v.*(Cy*u) + u*Dx' + Dy*u);
    v = v + dt*(u.*(v*Cx') + v.*(Cy*v) + v*Dx' + Dy*v);
    Xu = [Xu u(:)];
    Xv = [Xv v(:)];
end
%}
disp('FOM computed')
toc
%}

%% Import Data %%
%%{
load('snapshotmatrices_regularized_Re3000.mat')
X = [XU; XV];
%}
Xmean = mean(X,2);

%Xu = X(1:Nx*Ny,:);
%Xv = X((Nx*Ny+1):(2*Nx*Ny),:);
M_size = 2*Nx*Ny;

%save('FOM_data2D_100','Cx','Cy','Dx','Dy','D','CX','CY','Nx','Ny','Nt','dx','dt','timespan')

figure
for tt = 1:0.05*Nt:(Nt+1)
    %%{
    p = quiver(reshape(XU(:,tt),[Ny,Nx]),reshape(XV(:,tt),[Ny,Nx]),7); axis equal tight; set(gca,'Ydir','reverse');
    drawnow
    %{
    basename1 = 'LDC_U';
    filename1 = [basename1,num2str(tt),'.jpg'];
    saveas(p,filename1);
    %}
end
clear tt

for tt = 1:0.05*Nt:(Nt+1)
    %{
    p = pcolor(reshape(Xv(:,tt),[Ny,Nx])); colorbar; shading interp; colormap('jet'); axis equal; set(gca,'Ydir','reverse');
    
    basename1 = '2D_Burgers_V';
    filename1 = [basename1,num2str(tt),'.jpg'];
    saveas(p,filename1);
    %}
end
clear tt

%% POD Part %%
%X(:,1) = [];          % Special case: empty vector at beginning of simulation to be deleted
tr = 1.0;                       % Truncation time factor - trucates original Solution Matrix X
% y0 = X(:,1);          % Add or remove %-sign in, to suppress or enable subtraction of IC
y0 = mean(X,2);
% y0 = zeros(size(X(:,1))); 
skip = 1:2:tr*Nt+1; Nskip = length(skip);
Xtr = X(:,skip)-repmat(y0,1,Nskip);        % Snapshot matrix up to time tr*Nt, note: first column is initial condition
tic
[U,S,V] = svd(Xtr,'econ'); 
disp('SVD computed')
toc
[N,K] = size(Xtr);
k = 40;                         % Number of chosen POD modes

Uk = U(1:N,1:k);
Sk = S(1:k,1:k);
Vk = V(1:K,1:k);

Xlra = Uk*Sk*Vk'+repmat(y0,1,Nskip);               % Low-rank-approximation of Xtr
Vlra = Sk*Vk';                  % temporal modes over snapshot time (pseudo time)
frob = norm(Xtr+repmat(y0,1,Nskip)-Xlra,'fro')
POD_energy = (trace(S(1:min(N,K),1:min(N,K)))-trace(Sk))/trace(S(1:min(N,K),1:min(N,K)))

aPOD = Uk'*(X-repmat(y0,1,size(X,2)));  %predicted aPOD modes (exceeding snapshot time)

% Plot low-rank-approximated snapshot matrix
%{
Xlrau = Xlra(1:Nx*Ny,:);
Xlrav = Xlra((Nx*Ny+1):(2*Nx*Ny),:);
for tt = 1:round(0.1*(K)):K
    
    p = quiver(reshape(Xlrau(:,tt),[Ny,Nx]),reshape(Xlrav(:,tt),[Ny,Nx]),7); colorbar; shading interp; colormap('jet'); axis equal; set(gca,'Ydir','reverse');
    drawnow
    %{
    basename1 = 'ProjErrVec';
    filename1 = [basename1,num2str(tt),'.jpg'];
    saveas(p,filename1);
    %}
end
clear tt
%}

%%{
figure
plot(1:K,Vlra(1,:),1:K,Vlra(2,:),1:K,Vlra(3,:),1:K,Vlra(4,:),1:K,Vlra(5,:))
title('aPOD Modes up to cut off time')
legend('1st','2nd','3rd','4th','5th')
%}

%{
for pp = 1:k
    p = pcolor(transpose(reshape(Uk(:,pp),[Nx,Ny]))); axis equal; shading interp;
    
    %basename1 = '2D_ConvDiff_PODmode';
    %filename1 = [basename1,num2str(pp),'.jpg'];
    %saveas(p,filename1);
end
clear pp
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
QX = Uk' *qX;
QY = Uk' *qY;
Qrom = QX + QY;            % Completely precomputed quadratic matrix for nonlinear terms
clear ii jj qX qY QX QY

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
LX2 = Uk'* lX2;
LY2 = Uk'* lY2;
L2 = LX2 + LY2;
clear ii lX2 lY2 LX2 LY2

% LX3 & LY3
lX3 = zeros(M_size,k);
lY3 = zeros(M_size,k);
for ii = 1:k
    lX3(:,ii) = uu0(:,1).* CXUk(:,ii);
    lY3(:,ii) = vv0(:,1).* CYUk(:,ii);
end
LX3 = Uk'* lX3;
LY3 = Uk'* lY3;
L3 = LX3 + LY3;
clear ii lX3 lY3 LX3 LY3

% LX4 & LY4 
LX4 = Uk'* (uu0.* CXuv);
LY4 = Uk'* (vv0.* CYuv);
Crom = LX4 + LY4 + D0;
clear LX4 LY4 D0

Lrom = L2 + L3 + Dstar;          % Completely precomputed linear matrix for linear terms 
clear L2 L3 Dstar

%% Solving the ROM
t_max = 5.0;                              % total simulation time
dt = 0.001;                          % integration time size (not actual time discretization step size)
Na_t = t_max/dt;                      % number of temporal steps/ snapshots! 
timespan = 0:dt:t_max;             % integration timespan in [seconds]     
a0_ROM = Uk'*(X(:,1)-y0);                % initial conditions for a(t)

tic
%%{
[~,a_dot] = ode45(@(t,a) ode2DBurgersROM(a,Qrom,Lrom,Crom,Q3D,k), timespan, a0_ROM);
a_t = a_dot'; clear a_dot          % solution matrix for the temporal modes a(t) at different times
%}
%{
a_t = zeros(k,length(timespan));
a_t(:,1) = a0_ROM;
for tt = 1:Na_t
    a_t(:,tt+1) = a_t(:,tt) + dt*(Qrom * kron(a_t(:,tt),a_t(:,tt)) + Lrom*a_t(:,tt) + Crom); 
end
%}
disp('ROM computed')
toc

f_rom = Uk*a_t + y0;
u_rom = f_rom(1:Nx*Ny,:);
v_rom = f_rom((Nx*Ny+1):(2*Nx*Ny),:);

%save('VelocityROM_U_V','u_rom','v_rom');
%save('temporal_modes','aPOD','a_t');

%% Plotting/Visualization Section
%load('VelocityROM_U_V.mat')
figure
for tt = 1:round(0.05*Na_t):Na_t+1
    %%{
    p = quiver(reshape(u_rom(:,tt),[Ny,Nx]),reshape(v_rom(:,tt),[Ny,Nx]),7); colorbar; shading interp; colormap('jet'); axis equal; set(gca,'Ydir','reverse');
    drawnow
    %{
    basename1 = 'LDC_ROM_vec';
    filename1 = [basename1,num2str(tt),'.jpg'];
    saveas(p,filename1);
    %}
end
clear tt

%{
Dx_int = 1/nu*Dx(2:Nx-1,2:Nx-1);
Dy_int = 1/nu*Dy(2:Ny-1,2:Ny-1);
tic
tmp = full(kron(Dx_int,eye(Ny-2)) + kron(eye(Nx-2),Dy_int));
D_psi_int = (tmp)\eye((Nx-2)*(Ny-2));
disp('inverse laplacian computed')      %needs 218 seconds for a 128x128 grid
toc
clear tmp
Psi_rom = zeros(Ny,Nx);
Psi_rom_int = Psi_rom(2:Ny-1,2:Nx-1);
%}

%{
for tt = 1:round(0.05*Na_t):Na_t+1
    U_rom = reshape(u_rom(:,tt),[Ny,Nx]);
    V_rom = reshape(v_rom(:,tt),[Ny,Nx]);
    W_rom = -V_rom*Cx'+ Cy*U_rom;
    W_rom_int = W_rom(2:Ny-1,2:Nx-1);
    Psi_rom_int(:) = -D_psi_int*W_rom_int(:);   % compute psi on interior domain only through in inverse laplacian
    Psi_rom(2:Ny-1,2:Nx-1) = Psi_rom_int;
    
    p = pcolor(Psi_rom); colorbar; shading interp; colormap('jet'); axis equal; set(gca,'Ydir','reverse');
    drawnow
    %{
    basename1 = 'Re30000_ROM_k20_psi';
    filename1 = [basename1,num2str(tt),'.png'];
    saveas(p,filename1);
    %}
end
clear tt
%}

tt = (1:size(aPOD,2))/size(aPOD,2); ti = (1:size(a_t,2))/size(a_t,2);
figure
plot(tt,aPOD(1,:),'b-.',tt,aPOD(2,:),'r-.',tt,aPOD(3,:),'g-.',tt,aPOD(4,:),'c-.',tt,aPOD(5,:),'y-.',...
    ti,a_t(1,:),'b',ti,a_t(2,:),'r',ti,a_t(3,:),'g',ti,a_t(4,:),'c',ti,a_t(5,:),'y')
title('a_{ROM} vs a_{POD} LDC'); xlabel('time'); ylabel('a(t)'); 
legend('1st POD','2nd POD','3rd POD','4th POD','5th POD','1st ROM','2nd ROM','3rd ROM','4th ROM','5th ROM')