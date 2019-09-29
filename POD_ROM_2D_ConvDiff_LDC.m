%% Galerkin-POD ROM for 2D Convection-Diffusion in a time-averaged LDC flow field %%
clear
close all

%% Pre Setup %%
nu = 0.001;     % Scalar Diffusivity      
c_x = 1;        % Transport speed in x direction
c_y = 1;        % Transport speed in y direction   

L_x = 1;       % Domain Length in x-Direction
L_y = 1;        % Domain Length in y-Direction
Nx = 100;       
Ny = 100;
%M_size = Nx*Ny;
dx = L_x/(Nx-1) 
dy = L_y/(Ny-1)

tmax = 40;      % end time
Nt = 2000;       % number of timesteps
dt = tmax/Nt;      % timestep size 

%% Load flow field from Re=30000 LDC (and interpolate from 128x128 grid to current NX x NY grid)
load('LDC_100snapshots_Nx128_Re30000.mat')
disp('Turbulent flow field loaded')

[xq,yq] = meshgrid(0:dx:L_x);
[xold,yold] = meshgrid(0:1/127:1);
Ufield = interp2(xold,yold,reshape(mean(XU_100,2),[128,128]),xq,yq,'spline');
Vfield = interp2(xold,yold,reshape(mean(XV_100,2),[128,128]),xq,yq,'spline');
figure
quiver(Ufield,Vfield,3)
axis equal; set(gca,'Ydir','reverse');
title('Fig. 3: Time-averaged flow field, Vector Plot')

clear XU_100 XV_100 xold yold
Ufield = Ufield(2:Ny-1,2:Nx-1);
Vfield = Vfield(2:Ny-1,2:Nx-1);

%% Creation of system matrices for FOM %%

% ===== Boundary Conditions ===== 
%{
eastBC = ones(Ny-2,1);             % for Dirichlet BC's
westBC = ones(Ny-2,1);
northBC = ones(Nx-2,1);
southBC = ones(Nx-2,1);
%}
%%{
eastBC = zeros(Ny-2,1);             % for Dirichlet BC's
westBC = zeros(Ny-2,1);
northBC = zeros(Nx-2,1);
southBC = zeros(Nx-2,1);
%
% ===== Boundary Conditions =====

Dx = nu/dx^2*(diag(-2*ones(Nx,1)) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1));                  % x second derivative Matrix
Dx = sparse(Dx);
Dy = nu/dy^2*(diag(-2*ones(Ny,1)) + diag(ones(Ny-1,1),1) + diag(ones(Ny-1,1),-1));  % y second derivative Matrix
Dy = sparse(Dy);

Cx = c_x/(2*dx)*(diag(-1*ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1));                                       % x first derivative Matrix
Cx = sparse(Cx); 
Cy = c_y/(2*dy)*(diag(-1*ones(Ny-1,1),1) + diag(ones(Ny-1,1),-1));                       % y first derivative Matrix
Cy = sparse(Cy);

Dx_int = Dx(2:Nx-1,2:Nx-1);
Dy_int = Dy(2:Ny-1,2:Ny-1);
BDxEW = nu/(dx^2)*[westBC zeros(Ny-2,Nx-4) eastBC];
BDyNS = nu/(dy^2)*[northBC'; zeros(Ny-4,Nx-2); southBC'];

Cx_int = Cx(2:Nx-1,2:Nx-1);
Cy_int = Cy(2:Ny-1,2:Ny-1);
BCxEW = c_x/(2*dx)*[westBC zeros(Ny-2,Nx-4) -eastBC];
BCyNS = c_y/(2*dy)*[northBC'; zeros(Ny-4,Nx-2); -southBC'];

%% Define spatial functions %%
r1 = zeros(Ny,Nx);                   % spatial oscillation shape functions
r2 = zeros(Ny,Nx);
for jj = 1:Ny      
    for ii = 1:Nx
        xi = ii/Nx*L_x;
        yj = jj/Ny*L_y;
        
        r1(jj,ii) =  30*exp(-(((xi-0.80)/0.07)^8 + ((yj-0.20)/0.07)^8))+...
                    - 30*exp(-(((xi-0.75)/0.07)^8 + ((yj-0.25)/0.07)^8));
        r2(jj,ii) = 30*exp(-(((xi-0.80)/0.05)^8 + ((yj-0.80)/0.05)^8))+...
                      - 30*exp(-(((xi-0.85)/0.05)^8 + ((yj-0.85)/0.05)^8));

    end
end
clear ii jj

surf(r1 + r2); axis equal; set(gca,'Ydir','reverse')
title('Fig. 1: Surface plot of spatial functions')
r1 = r1(2:Ny-1,2:Nx-1);
r2 = r2(2:Ny-1,2:Nx-1);

f1 = 1;                 % Frequency of harmonic excitation for r1
f2 = 2;                 % Frequency of harmonic excitation for r2

% Initialize Noise
%{
noise = zeros((Nx-2)*(Ny-2),Nt+1);
for tt = 1:Nt+1
     tmp = (r1+r2).*(20*(rand(size(r1))-0.5)) ;
     noise(:,tt) = tmp(:);
end
%}
load('noise2.mat')
figure
contourf(reshape(noise(:,1),[Ny-2,Nx-2])); colorbar; axis equal; set(gca,'Ydir','reverse');
title('Fig. 2: Position of spatial functions')
disp('New noise initialized')

%% Solve the FOM
timespan = 0:dt:tmax;             % integration timespan in [seconds]
Asys = 10*Ufield(:).*kron(Cx_int,eye(Nx-2)) + 10*Vfield(:).*kron(eye(Ny-2),Cy_int) + kron(Dx_int,eye(Nx-2)) + kron(eye(Ny-2),Dy_int);
bsys  = 10*Ufield(:).*BCxEW(:) + 10*Vfield(:).*BCyNS(:) + BDxEW(:) + BDyNS(:);

IC_FOM = max(northBC)*ones((Ny-2)*(Nx-2),1);
tFOM = tic;
[~,y] = ode45(@(t,y) ode2DConvDiffMatrix(t,y,Asys,bsys,Nx,Ny,r1(:),r2(:),f1,f2,dt,noise),timespan,IC_FOM);
Xnoise = y'; clear y                                % solution matrix ==> for snapshot matrix
tFOM = toc(tFOM)

figure
for tt = 0.1*Nt:10:(Nt+1)
    pcolor(reshape(Xnoise(:,tt),[Ny-2,Nx-2])); axis equal; set(gca,'Ydir','reverse'); shading interp; colorbar
    title(['Conv-Diff FOM at timestep = ' num2str(tt)]); caxis([-2 4]);
    drawnow 
    %{
    basename1 = '2D_ConvDiff';
    filename1 = [basename1,num2str(tt),'.jpg'];
    saveas(p,filename1);
    %}
end
clear tt

%{
figure % This draws 9 instances of the transient
ctr = 1;
for tt = [40:40:360]
    subplot(3,3,ctr); ctr=ctr+1;
    pcolor(reshape(Xnoise(:,tt),[Ny-2,Ny-2])); axis equal; set(gca,'Ydir','reverse'); shading interp; colorbar
    title(['Conv-Diff FOM at timestep = ' num2str(tt)]); caxis([-2 4]);
    drawnow 
end
clear tt
%}

%% END OF FULL-ORDER MODEL (FOM)
close all
clear BCxEW BCyNS BDxEW BDyNS Cx c_x c_y Cx_int Cy Cy_int Dx Dx_int Dy Dy_int northBC eastBC southBC westBC Ufield Vfield IC_FOM tmp

%% POD algorithm
tr = 1;                       % Truncation time factor - trucates original Solution Matrix X
y0 = zeros((Ny-2)*(Nx-2),1);          % Add or remove percent sign in order to suppress or enable subtraction of IC 
X0 = repmat(y0,1,length(1:(tr*Nt+1)));
Xtr = Xnoise(:,1:(tr*Nt+1))-X0;        % Snapshot matrix up to time tr*Nt, note: first column is initial condition
tic
[U,S,V] = svd(Xtr,'econ');  disp('SVD computed')
toc 
[N,K] = size(Xtr);
k = 200;                         % Number of chosen POD modes

Uk = U(1:N,1:k);
Sk = S(1:k,1:k);
Vk = V(1:K,1:k);

Xlra = Uk*Sk*Vk'+X0;               % Low-rank-approximation of Xtr
Vlra = Sk*Vk';                     % a_POD(t)
frob = norm(Xtr + X0 - Xlra,'fro')/norm(Xtr + X0,'fro')
unused_POD_energy = (trace(S(1:min(N,K),1:min(N,K))).^2-trace(Sk).^2)/trace(S(1:min(N,K),1:min(N,K))).^2

figure
for kk = 1:k
    loglog(kk,Sk(kk,kk),'*'), hold on
end
loglog(1:min(N,K),diag(S),'k:')
xlabel('Mode i'), ylabel('\sigma'), xlim([0 1000])
legend('\sigma_1','\sigma_2','\sigma_3','\sigma_4','\sigma_5','\sigma_6','\sigma_7','\sigma_8', ...
'\sigma_9','\sigma_{10}','\sigma_{11}','\sigma_{12}','\sigma_{13}','\sigma_{14}','\sigma_{15}'), set(gca,'FontSize',16), grid on

% Plot low-rank-approximated snapshot matrix
%{
figure
for tt = 0.1*Nt:10:(Nt+1)
    subplot(1,2,1)
    pcolor(reshape(Xnoise(:,tt),[Ny-2,Nx-2])); axis equal; set(gca,'Ydir','reverse'); shading interp; colorbar
    title(['Conv-Diff FOM at timestep = ' num2str(tt)]); caxis([-2 4]);    
    
    subplot(1,2,2)
    pcolor(reshape(Xlra(:,tt),[Ny-2,Nx-2])); axis equal; set(gca,'Ydir','reverse'); shading interp; colorbar
    title(['Low-rank k=' num2str(k) ' reconstruction at timestep = ' num2str(tt)]); caxis([-2 4]);    
    drawnow 
end; clear tt
%}

%% Galerkin Projection & ROM part
A_rom = Uk'*Asys*Uk;            % GALERKIN PROJECTION - obtain the kxk matrix to solve for a_t(1) to a_t(k)
eigval = eig(A_rom)             % Check Eigenvalues for potentially unstable behavior 
b_rom = Uk'*bsys;
c_rom = Uk'*Asys*y0;
r1_rom = Uk'*r1(:);
r2_rom = Uk'*r2(:);
noise_rom = Uk'*noise;
     
a0_rom = Uk'*Xnoise(:,1); 

tROM = zeros(1,10);
for tt = 1:10
    tmp = tic;
    [t,y] = ode45(@(t,y) ode2DConvDiffPODROM(t,y,A_rom,b_rom,c_rom,Nx,Ny,r1_rom,r2_rom,f1,f2,dt,noise_rom), timespan, a0_rom);
    tROM(tt) = toc(tmp);
    %disp('POD-ROM computed')
end
tROM_avg = mean(tROM)
a_rom = y'; clear y  % solution matrix for the temporal modes a(t) at different times

X_rom = Uk*a_rom + X0(:,1:length(a_rom));

% Plot evolution of a_POD vs a_ROM over time
tt = timespan;
figure
plot(tt,Vlra(1,:),'b-.',tt,Vlra(2,:),'r-.',tt,Vlra(3,:),'g-.',tt,Vlra(4,:),'m-.', ...
    tt,a_rom(1,:),'b',tt,a_rom(2,:),'r',tt,a_rom(3,:),'g',tt,a_rom(4,:),'m')
%title('a_{POD} vs a_{ROM} over time'), 
xlim([0,10])
legend('1st a_{POD}','2nd a_{POD}','3rd a_{POD}','4th a_{POD}',...
'1st a_{ROM}','2nd a_{ROM}','3rd a_{ROM}','4th a_{ROM}')
xlabel('t'), ylabel('a(t)'), grid on, set(gca,'FontSize',16)

TKE_POD = 0.5*sum((U'*Xtr).^2); tt = timespan;
figure
plot(tt,TKE_POD), xlabel('t'), ylabel('TKE'), grid on, set(gca,'FontSize',16)
%title('TKE_{POD} over time')

% Plot ROM vs FOM
figure
for tt = 0.1*Nt:10:(Nt+1)
    subplot(1,2,1)
    pcolor(reshape(Xnoise(:,tt),[Ny-2,Nx-2])); axis equal; set(gca,'Ydir','reverse'); shading interp; colorbar
    title(['Conv-Diff FOM at timestep = ' num2str(tt)]); caxis([-3 3]);    
    
    subplot(1,2,2)
    pcolor(reshape(X_rom(:,tt),[Ny-2,Nx-2])); axis equal; set(gca,'Ydir','reverse'); shading interp; colorbar
    title(['ROM k=' num2str(k) ' reproduction at timestep = ' num2str(tt)]); caxis([-3 3]);    
    drawnow 
end
clear tt

%% Compute Errors
%{
rel_err_inf = norm(Xnoise-X_rom,'inf')/norm(Xnoise,'inf')
rel_err_fro = norm(Xnoise-X_rom,'fro')/norm(Xnoise,'fro')
rel_err_2 = norm(Xnoise-X_rom,2)/norm(Xnoise,2)
%}
%%{
rel_err_inf = norm(Xnoise(:,0.1*Nt:end)-X_rom(:,0.1*Nt:end),'inf')/norm(Xnoise(:,0.1*Nt:end),'inf')
rel_err_fro = norm(Xnoise(:,0.1*Nt:end)-X_rom(:,0.1*Nt:end),'fro')/norm(Xnoise(:,0.1*Nt:end),'fro')
%rel_err_2 = norm(Xnoise(:,0.1*Nt:end)-X_rom(:,0.1*Nt:end),2)/norm(Xnoise(:,0.1*Nt:end),2)
%}
rel_err_2 = zeros(1,size(Xnoise,2));
for tt = 0.1*Nt:length(rel_err_2)
    rel_err_2(tt) = norm(Xnoise(:,tt)-X_rom(:,tt))/norm(Xnoise(:,tt));
end 
mean(rel_err_2,2)

%% Create film
writer_obj = VideoWriter('2D_SPOD_ROM.avi');
writer_obj.FrameRate = 20;
open(writer_obj);
% write the frames to the video
for i=1:length(Film)
    % convert the image to a frame
    frame = Film(i) ;    
    writeVideo(writer_obj, frame);
end
% close the writer object
close(writer_obj);
