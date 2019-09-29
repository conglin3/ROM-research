%% Lid Driven Cavity in Vorticity Stream Function Form %%
clear
close all

%% Pre Setup %%
nu = 1.0/100;             

L_x = 1;       % Domain Length in x-Direction
L_y = 1;        % Domain Length in y-Direction
Nx = 50;      
Ny = 50;
dx = L_x/(Nx-1)
dy = L_y/(Ny-1)

%% Creation of system matrices for FOM %%

Dx = nu*1/(dx^2)*(diag(-2*ones(Nx,1)) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1));                  % x second derivative Matrix
Dx = sparse(Dx);
Dy = nu*1/(dy^2)*(diag(-2*ones(Ny,1)) + diag(ones(Ny-1,1),1) + diag(ones(Ny-1,1),-1));  % y second derivative Matrix
Dy = sparse(Dy);

Cx = 1/(2*dx)*(diag(ones(Nx-1,1),1) + diag(-1*ones(Nx-1,1),-1));                                       % x first derivative Matrix
Cx = sparse(Cx); 
Cy = 1/(2*dy)*(diag(ones(Ny-1,1),1) + diag(-1*ones(Ny-1,1),-1));                       % y first derivative Matrix
Cy = sparse(Cy);

%Bns = zeros(Ny,Ny); Bns(1,1) = 1; Bns(1,2) = -1; Bns(Ny,Ny) = 1; Bns(Ny,Ny-1) = -1; Bns = 2/(dy^2)*Bns;
%Bew = zeros(Nx,Nx); Bew(1,1) = 1; Bew(2,1) = -1; Bew(Nx,Nx) = 1; Bew(Nx-1,Nx) = -1; Bew = 2/(dx^2)*Bew;
%Bu = zeros(Ny,Nx); Bu(1,:) = ones(1,Nx);

Dx(1,:) = nu*1/(dx^2)*[[2 -5 4 -1] zeros(1,Nx-4)];
Dx(Nx,:)= nu*1/(dx^2)*[zeros(1,Nx-4) [-1 4 -5 2]];
Dy(1,:) = nu*1/(dy^2)*[[2 -5 4 -1] zeros(1,Ny-4)];
Dy(Ny,:)= nu*1/(dy^2)*[zeros(1,Ny-4) [-1 4 -5 2]];
Cx(1,:) = 1/(2*dx)* [-3 4 -1 zeros(1,Nx-3)];
Cx(Nx,:) = 1/(2*dx)* [zeros(1,Nx-3) 1 -4 3];
Cy(1,:) = 1/(2*dy)* [-3 4 -1 zeros(1,Ny-3)];
Cy(Ny,:) = 1/(2*dy)* [zeros(1,Ny-3) 1 -4 3];

Dx_int = 1/nu*Dx(2:Nx-1,2:Nx-1);
Dy_int = 1/nu*Dy(2:Ny-1,2:Ny-1);
tic
tmp = full(kron(Dx_int,eye(Ny-2)) + kron(eye(Nx-2),Dy_int));
D_psi_int = (tmp)\eye((Nx-2)*(Ny-2));
disp('inverse laplacian computed')      %needs around 120 seconds for a 128x128 grid
toc
clear tmp

%Cx_int = Cx(2:Nx-1,2:Nx-1);
%Cy_int = Cy(2:Ny-1,2:Ny-1);

%% Solving the FOM %%
%%{
w = zeros(Ny,Nx);             
psi = zeros(Ny,Nx);            % initial fields all zero
U = zeros(Ny,Nx);
V = zeros(Ny,Nx);
%}

t_max = 10;                              % total simulation time
dt_int = 0.01;                         % time integration step size
N_t = t_max/dt_int;                     % number of temporal steps/ snapshots! 
timespan = 0:dt_int:t_max;              % integration timespan in [seconds]

fsnap = 1;
Nsnap = N_t/fsnap; 
tsnap = 1;

Xpsi = zeros(Nx*Ny,Nsnap+1); Xpsi(:,1) = psi(:);
Xw = zeros(Nx*Ny,Nsnap+1);   Xw(:,1) = w(:);
XU = zeros(Ny*Nx,Nsnap+1);   XU(:,1) = U(:);
XV = zeros(Ny*Nx,Nsnap+1);   XV(:,1) = V(:);

w_int = w(2:Ny-1,2:Nx-1);
psi_int = psi(2:Ny-1,2:Nx-1);

U_north = 2.0*exp(-((( dx*[0 meshgrid(1:Nx-1,1)]-0.5)/0.4).^10));
%U_north = -(dx*[0 meshgrid(1:Nx-1,1)]-1).^2 + 1;
%U_north = -1*ones(1,Nx);
%disp(['Re = ',num2str(U_north(1,Nx/2)*L_x/nu)])
%First run: sin(0.15*pi*timespan(tt).^2) with U_north=10 and nu = 0.001
%Second run: sin(0.01*pi*timespan(tt).^3 + pi/4)
%Third run: -cos(0.0015*pi*(timespan+t_max).^2 + 0.6*pi)

tic
for tt = 1:N_t+1
    if(mod(tt,100) == 0) 
        if(max(max(w_int))> 1e6)
            disp('diverging')
            break
        end
        disp(['timestep ',num2str(tt)])
        pcolor(w); colorbar; shading interp; colormap('jet'); axis equal tight;
        xlabel('x points');ylabel('y points');caxis([-10 10]); set(gca,'Ydir','reverse'); 
        title(['Vorticity contours, Re=' num2str(max(U_north)/nu) ' timestep=' num2str(tt)]);drawnow
    end    
    
    psi_int(:) = -D_psi_int*w_int(:);   % compute psi on interior domain only through inverse laplacian
    psi(2:Ny-1,2:Nx-1) = psi_int;       
%%{   
    w(1,:) = (-psi(2,:))*2/(dy^2) + U_north*2/dy;       % North BC psi(1,:)=0       
    w(:,1) = (-psi(:,2))*2/(dx^2);                      % West BC psi(:,1)=0
    w(Ny,:) = (-psi(Ny-1,:))*2/(dy^2);                 % South BC psi(Ny,:)=0
    w(:,Nx) = (-psi(:,Nx-1))*2/(dx^2);                 % East BC psi(:,Nx)=0
%}     
%{
    w(1,:) = (psi(4,:)-11*psi(2,:))/(dy^2) + sin(0.01*pi*timespan(tt).^3 + pi/4)*U_north*8/dy;       % North BC psi(1,:)=0
    w(:,1) = (psi(:,4)-11*psi(:,2))/(dx^2);                      % West BC psi(:,1)=0
    w(Ny,:) = (psi(Ny-3,:)-11*psi(Ny-1,:))/(dy^2);                 % South BC psi(Ny,:)=0
    w(:,Nx) = (psi(:,Nx-3)-11*psi(:,Nx-1))/(dx^2);                 % East BC psi(:,Nx)=0
%}    
    U = Cy*psi;
    U(1,:) = U_north; 
    U(Ny,:) = zeros(1,Nx);
    V = -psi*Cx';
    V(:,1) = zeros(Ny,1);
    V(:,Nx) = zeros(Ny,1);
    
    tmp = (-U).*(w*Cx') + (-V).*(Cy*w) + (w*Dx' + Dy*w);          % temporary RHS with incorrect BCs but correct interior
    w_int = w(2:Ny-1,2:Nx-1) + dt_int*(tmp(2:Ny-1,2:Nx-1));      % updated vorticity interior after one time step integration
    w(2:Ny-1,2:Nx-1) = w_int;
    
%%{
    if(mod(tt,fsnap)==0)
        Xw(:,tsnap+1) = w(:);          % Vorticity snapshots
        Xpsi(:,tsnap+1) = psi(:);      % Stream function snapshots
        XU(:,tsnap+1) = U(:);          % U velocity snapshots
        XV(:,tsnap+1) = V(:);          % V velocity snapshots
        tsnap = tsnap+1;
    end
%}
%{
    if(mod(tt,4)==0)
        Xw = [Xw w(:)];          % Vorticity snapshots
        Xpsi = [Xpsi psi(:)];      % Stream function snapshots
        XU = [XU U(:)];          % U velocity snapshots
        XV = [XV V(:)];          % V velocity snapshots 
    end
%}
end
toc
clear tt

%% Save and visualize
%%{
Xpsi(:,1) = [];
Xw(:,1) = [];
XU(:,1) = [];
XV(:,1) = [];
%}
save('snapshotmatrices_Re200_10sec','Xw','Xpsi','XU','XV');

figure
for tt = 1:round(0.005*size(Xw,2)):size(Xw,2)
    %%{
    p = pcolor(reshape(Xw(:,tt),[Ny,Nx])); colorbar; shading interp; colormap('jet'); axis equal tight;
    xlabel('x points');ylabel('y points');caxis([-8 8]); set(gca,'Ydir','reverse'); 
    title(['Vorticity contours, Re=' num2str(max(U_north)/nu) ' timestep=' num2str(tt)]);drawnow
    %{
    basename1 = 'LDC_w';
    filename1 = [basename1,num2str(tt),'.png'];
    saveas(p,filename1);
    %}
end
clear tt

for tt = 1:round(0.005*size(XU,2)):size(XU,2)
    %%{
    p = quiver(reshape(XU(:,tt),[Ny,Nx]),reshape(XV(:,tt),[Ny,Nx]),7); axis equal tight; set(gca,'Ydir','reverse');
    drawnow
    %}
    %{
    basename1 = 'LDC_vec';
    filename1 = [basename1,num2str(tt),'.png'];
    saveas(p,filename1);
    %}
end
clear tt

for tt = 1:round(0.005*size(Xpsi,2)):size(Xpsi,2)
    %%{
    p = pcolor(reshape(Xpsi(:,tt),[Ny,Nx])); colorbar; shading interp; colormap('jet'); axis equal;
    xlabel('x points');ylabel('y points');set(gca,'Ydir','reverse'); %caxis([-0.01 0.05]);
    title(['Streamline contours, Re=' num2str(max(U_north)/nu) ' timestep=' num2str(tt)]);drawnow
    %{
    basename1 = 'LDC_psi';
    filename1 = [basename1,num2str(tt),'.png'];
    saveas(p,filename1);
    %} 
end
clear tt

%{
% SOR iterative algorithm for the streamfunction %

MaxIt = 100;
MaxErr = 0.001;
beta = 1.5;                             %relaxation factor
gamma = 1/(2/dx^2 + 2/dy^2);

    for it = 1:MaxIt
        tmp = psi;        
        for ii = 2:Nx-1
            for jj = 2:Ny-1                 
                psi(jj,ii) = beta*gamma*( 1/(dy^2)*(psi(jj+1,ii) + psi(jj-1,ii)) + 1/(dx^2)*(psi(jj,ii+1) + psi(jj,ii-1)) + w(jj,ii)) + (1-beta)*psi(jj,ii);
            end
        end        
        err = 0;
        err = sumabs(tmp-psi);
        if err <= MaxErr 
            break 
        end    
    end
    clear ii jj
%}

