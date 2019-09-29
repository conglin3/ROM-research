%% SPOD-ROM for 2D Convection-Diffusion in a time-averaged LDC flow field %%
clear
close all
cd C:\Users\clin\Documents\MATLAB\_Codes_and_Data\2D_SPOD_ROM

%% Pre Setup %%
nu = 0.001;     % Scalar Diffusivity      
c_x = 1;        % Transport speed for Cx matrix
c_y = 1;        % Transport speed for Cy matrix   

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

clear XU_100 XV_100 xold yold
Ufield = Ufield(2:Ny-1,2:Nx-1); xq = xq(2:Ny-1,2:Nx-1); 
Vfield = Vfield(2:Ny-1,2:Nx-1); yq = yq(2:Ny-1,2:Nx-1); 

figure
%quiver(Ufield,Vfield,7)
pcolor(xq,yq,((10*Vfield).^2+(10*Ufield).^2).^(0.5))
axis equal tight; set(gca,'Ydir','reverse'); shading interp; caxis([0 3.5]), colorbar
%title('Fig. 3: Time-averaged flow field, Vector Plot')
set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dy L_y-dy],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01],'XAxisLocation','top','FontSize',16)

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
%}
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

%% Initialize the Spatial functions %%
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
clear ii jj xi yj

surf(1:Nx,1:Ny,(r1+r2)); axis equal; set(gca,'Ydir','reverse');colorbar;
%title('Fig. 1: Surface plot of spatial functions')
set(gca,'Ydir','reverse','XTick',[1 Nx],'XTickLabel',[0 1],'YTick',[1 Ny],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01],'XAxisLocation','top','FontSize',16)
r1 = r1(2:Ny-1,2:Nx-1);
r2 = r2(2:Ny-1,2:Nx-1);

f1 = 1;                 % Frequency of harmonic excitation for r1
f2 = 2;                 % Frequency of harmonic excitation for r2
%{
noise = zeros((Nx-2)*(Ny-2),Nt+1);
for tt = 1:Nt+1
     tmp = (r1+r2).*(20*(rand(size(r1))-0.5)) ;
     noise(:,tt) = tmp(:);
end; disp('New noise initialized')
%}
load('noise2.mat'); disp('Noise loaded')
figure
contourf(xq,yq,reshape(noise(:,1),[Ny-2,Nx-2]),16), colorbar; axis equal tight; set(gca,'Ydir','reverse');
%title('Fig. 2: Spatial functions overlaid with one instance of noise')
set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dy L_y-dy],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01],'XAxisLocation','top','FontSize',16)
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
for tt = 0.1*Nt:1:(Nt+1)
    pcolor(reshape(Xnoise(:,tt),[Ny-2,Nx-2])); axis equal tight; set(gca,'Ydir','reverse'); shading interp; colorbar
    title(['Conv-Diff FOM at timestep = ' num2str(tt)]); caxis([-3 3]); xlabel('x');ylabel('y');
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
for tt = [11:20:171]
    subplot(3,3,ctr); ctr=ctr+1;
    pcolor(xq,yq,reshape(Xnoise(:,tt),[Ny-2,Ny-2])); axis equal tight; shading interp; colorbar
    title(['HFM at sec = ' num2str((tt-1)*dt)]); caxis([-3 3]); %xlabel('x');ylabel('y'); 
    set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dx L_x-dx],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01],'FontSize',12)%,'XAxisLocation','top')
    drawnow 
end
clear tt
%}

%% END OF FULL-ORDER MODEL (FOM)
close all 
clear BCxEW BCyNS BDxEW BDyNS Cx c_x c_y Cx_int Cy Cy_int Dx Dx_int Dy Dy_int northBC eastBC southBC westBC Ufield Vfield IC_FOM tmp
clear Q0 Q0hat_fk Q0hat Q_rom Q_rom_real Qhat_rom Qhat_fk Qhat_rom_fk Fhat_fk Film

%% SPOD algorithm
X = Xnoise(:,round(0.1*Nt):end);   % choose how many snapshots to be used for SPOD
%X = X - mean(X,2);          % optional: subtract mean --> THIS is what the original SPOD code by Schmidt et al. does
Nt = size(X,2);             % update number of timesteps actually used
Nx = size(X,1);             % Number of spatial DoFs for SPOD algorithm 

Nf = 1000;                  % Number of frequencies in a block
No = 942;                   % Block overlap number
Nb = floor((Nt-No)/(Nf-No)) % effective number of blocks

%Partition the snapshot matrix into Nb blocks 
Q = zeros(Nx,Nf,Nb); %Q0 = zeros(Nx,Nb);
for n = 1:Nb 
    Q(:,:,n) = X(:,(1+(n-1)*(Nf-No)):(Nf+(n-1)*(Nf-No)));
    %Q0(:,n) = mean(Q(:,:,n),2);
end
disp(['Matrix partitioned into Nb = ' num2str(Nb) ' blocks'])

%Define Hamming Window to reduce spectral leakage due to non-periodicity of the data
%window = 0.54-0.46*cos(2*pi*(0:Nf-1)/(Nf-1)); disp('Using hamming window for FFT')
%Or use uniform window for FFT
window = ones(1,Nf);   disp('Using uniform window for FFT')
winWeight = 1/mean(window);

%FFT of the Data 
Qhat = zeros(Nx,Nf,Nb); Q0hat = zeros(Nx,Nf,Nb);
for n = 1:Nb
    Qhat(:,:,n) = winWeight/Nf*fft(squeeze(Q(:,:,n)).*window,Nf,2);
    %Q0hat(:,:,n) = winWeight/Nf*fft(repmat(Q0(:,n),1,Nf).*window,Nf,2);
    disp(['Block ' num2str(n) ' FFT computed'])
end

% Plot some of the realizations in the frequency domain 
%%{
Lt = Nf;                    % Length of time signal
t = (0:1:Lt)*dt;            % Time vector
f=Nf/2*linspace(-1/t(end),1/t(end),Nf+1);   %two sided frequency vector
f(end)=[];

%%{
checkmat = zeros(sqrt(Nx),sqrt(Nx));
for jj = [0.2 0.8]      
    for ii = [0.2 0.8]      
        checkmat(round(ii*sqrt(Nx)),round(jj*sqrt(Nx))) = 1;
    end
end; clear ii jj
checkvec = find(checkmat(:)'); checkvec = reshape(checkvec,[sqrt(length(checkvec)),sqrt(length(checkvec))])'; checkvec = checkvec(:)';
figure
pcolor(xq,yq,checkmat); axis equal tight, grid on,
set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dx L_x-dx],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01],'XAxisLocation','top','FontSize',16)
%}

figure
ctr = 1;
for ii = checkvec
    subplot(3,3,ctr); ctr = ctr+1;
    plot(f,fftshift(squeeze(abs(Qhat(ii,:,1)-Q0hat(ii,:,1)))))
    xlim([-20,20]), xlabel('f'), ylabel('Amplitude'), grid on
    title(['DoF #' num2str(ii) ' in the Frequency Domain'])
    drawnow
end
%}

%Group same frequency fourier modes from all realizations into Q_fk matrices 
f_half = Nf/2*linspace(0,1/t(end),Nf/2+1);
nfk = length(f_half)
Qhat_fk = permute(Qhat(:,1:nfk,:),[1 3 2]);
%Q0hat_fk = permute(Q0hat(:,1:nfk,:),[1 3 2]);
%SINCE INPUT DATA IS REAL: correct Fourier coefficients for one-sided spectrum such that amplitudes match up
%Qhat_fk(:,:,2:end-1) = 2*Qhat_fk(:,:,2:end-1);
%Q0hat_fk(:,:,2:end-1) = 2*Q0hat_fk(:,:,2:end-1);
disp('Permuted and corrected for 1-sided spectrum')

%And POD on Qhat_fk for all the fk frequencies
Uspod = zeros(Nx,Nb,nfk);
Aspod = zeros(Nb,Nb,nfk);
Lspod = zeros(Nb,nfk);
tic
for fk = 1:nfk
    %{
    % This is Schmidt's original algorithm
    M_fk = squeeze(Qhat_fk(:,:,fk))'*squeeze(Qhat_fk(:,:,fk))/Nb;   %normalize by factor kappa
    [Theta_fk,Lambda_fk] = eig(M_fk);                               %careful! eig() does not order eigenvalues in any particular order
    
    [lambda,idx]        = sort(diag(Lambda_fk),'descend');          %hence sort eigenvalues 
    Theta_fk            = Theta_fk(:,idx);                          %...and columns of Theta matrices
    
    Uspod(:,:,fk) = squeeze(Qhat_fk(:,:,fk))*Theta_fk*diag(1./sqrt(lambda))/sqrt(Nb);
    Aspod(:,:,fk) = diag(sqrt(lambda))*Theta_fk';
    Lspod(:,fk) = abs(lambda);
    disp(['SPOD computed for f_' num2str(fk) ' = ' num2str(f_half(fk))])
    %}
    %%{
    [U,S,V] = svd(squeeze((Qhat_fk(:,:,fk))/sqrt(Nb)),'econ');
    %[U,S,V] = svd(squeeze((Qhat_fk(:,:,fk)-Q0hat_fk(:,:,fk))/sqrt(Nb)));  
    disp(['SVD computed for f_{' num2str(fk) '} = ' num2str(f_half(fk))])
    Uspod(:,:,fk) = U(:,1:Nb);
    Aspod(:,:,fk) = S(1:Nb,:)*V';
    Lspod(:,fk) = diag(S.^2);
    %}
end
toc; clear U S V

figure
semilogy(f_half,Lspod'), grid on
xlabel('f'), ylabel('SPOD modal energy')
legend('\Lambda_1(f)','\Lambda_2(f)','\Lambda_3(f)','\Lambda_4(f)','\Lambda_5(f)','\Lambda_6(f)','\Lambda_7(f)','\Lambda_8(f)',...
'\Lambda_9(f)','\Lambda_{10}(f)','\Lambda_{11}(f)','\Lambda_{12}(f)','\Lambda_{13}(f)','\Lambda_{14}(f)','\Lambda_{15}(f)')
set(gca,'FontSize',16); xlim([0,15])
%title('SPOD modal energy spectra: harmonically+randomly forced passive scalar transport')

%evaluate decomposition error (infinity norm)
err_data=zeros(nfk,1);
for fk = 1:nfk
    err_data(fk,1) = norm((Qhat_fk(:,:,fk))/sqrt(Nb) - Uspod(:,:,fk)*Aspod(:,:,fk),'inf');
end 
disp(['Maximum decomposition error (infinity norm) is ' num2str(max(err_data))])

figure
count = 1;
for fi = [1 21 41]
    for mi = [1 2 3]
        subplot(3,3,count)
        pcolor(xq,yq,reshape(abs(Uspod(:,mi,fi)),[sqrt(Nx),sqrt(Nx)])); shading interp, axis equal tight, %caxis(max(abs(caxis))*[-1 1]) 
        title(['f_{' num2str(fi) '}=' num2str(f_half(fi),'%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(Lspod(mi,fi),'%.2g')])
        set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dx L_x-dx],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01])%,'XAxisLocation','top')
        drawnow 
        count = count + 1;
    end
end 

%% Galerkin Projection & ROM part
%{
noise2 = zeros((Nx-2)*(Ny-2),Nt+1); %Optional: create new noise
for tt = 1:Nt+1
     tmp = (r1+r2).*(20*(rand(size(r1))-0.5)) ;
     noise2(:,tt) = tmp(:);
end
%}
clear Q0 Q0hat_fk Q0hat Q_rom Q_rom_real Qhat_rom Qhat_fk Qhat_rom_fk Fhat_fk Film

tt = timespan(length(timespan)-Nt+1:end);
%force = noise(:,length(timespan)-Nt+1:end)+bsys;
%force = r1(:)*cos(2*pi*f1*tt)+r2(:)*cos(2*pi*f2*tt)+bsys;
force = bsys + r1(:)*cos(2*pi*f1*tt)+r2(:)*cos(2*pi*f2*tt)+noise(:,length(timespan)-Nt+1:end);

%Partition the forcing into Nb blocks and window-fft these blocks
F = zeros(Nx,Nf,Nb); Fhat = zeros(Nx,Nf,Nb);
for n = 1:Nb
    F(:,:,n) = force(:,(1+(n-1)*(Nf-No)):(Nf+(n-1)*(Nf-No)));
    Fhat(:,:,n) = winWeight/Nf*fft(squeeze(F(:,:,n)).*window,Nf,2);
    %Fhat(:,1:2,n) = 0; Fhat(:,end,n) = 0;
end
disp(['Forcing partitioned into Nb = ' num2str(Nb) ' blocks'])

%Group same frequency forcing into matrices
Fhat_fk = permute(Fhat(:,1:nfk,:),[1 3 2]);
%For one-sided spectrum: Double the amplitudes
%Fhat_fk(:,:,2:end-1) = 2*Fhat_fk(:,:,2:end-1);
 
bk = Nb
A_rom = zeros(bk,Nb,nfk);
Qhat_rom_fk = zeros(Nx,Nb,nfk);
iL_rom_fk = zeros(bk,bk,nfk);
F_rom_fk = zeros(bk,Nb,nfk);
tic
for fk = 1:nfk
    Psi_fk = squeeze(Uspod(:,1:bk,fk));
    
    L_rom = Psi_fk'*(Asys)*Psi_fk;
    iL_rom_fk(:,:,fk) = 1i*2*pi*f_half(fk)*eye(bk,bk)-L_rom;
    
    F_rom_fk(:,:,fk) = Psi_fk'*Fhat_fk(:,:,fk);
end
disp('ROM matrix precomputations finished')
toc

tROM = zeros(1,1000);
for tt = 1:1000
    A_rom = zeros(bk,Nb,nfk);
    tmp = tic;
    for fk = 1:nfk
        A_rom(:,:,fk) = 1/sqrt(Nb)*(squeeze(iL_rom_fk(:,:,fk))\squeeze(F_rom_fk(:,:,fk)));
        %disp(['SPOD-ROM computed for f_' num2str(fk)])
    end
    tROM(tt) = toc(tmp);
    %disp('ROM computed: A_rom evolution modes')
end
tROM_avg = mean(tROM)

tic
for fk = 1:nfk
    Qhat_rom_fk(:,:,fk) = sqrt(Nb)*squeeze(Uspod(:,1:bk,fk))*squeeze(A_rom(:,:,fk));    % Back-projection into original space
end
disp('Back-projection into N-dim space finished')
toc

clear F Fhat Psi_fk L_rom iL_rom_fk F_rom_fk 

%{
err_rom = zeros(nfk,1);
for fk = 1:nfk
    err_rom(fk,1) = norm(Aspod(1:kb,:,fk) - A_rom(:,:,fk))/norm(Aspod(1:kb,:,fk));
end
max(err_rom)
%}

%Convert back from single sided spectrum to double sided (for iDFT purposes)
%Qhat_rom_fk(:,:,2:end-1) = 0.5*Qhat_rom_fk(:,:,2:end-1);
Qhat_rom = permute(Qhat_rom_fk(:,:,:),[1 3 2]);
Qhat_rom = [Qhat_rom(:,:,:) conj(Qhat_rom(:,(end-1):-1:2,:))];
disp('Converted back to double sided DFT spectrum')

figure; ctr = 1;
for ii = checkvec
    subplot(2,2,ctr); ctr = ctr+1;
    semilogy(f,fftshift(squeeze(abs(Qhat(ii,:,1)))),f,fftshift(squeeze(abs(Qhat_rom(ii,:,1)))),'-.')
    legend('Data','SPOD-ROM'), grid on
    xlim([-15,15]), xlabel('f'), %ylabel('Amplitude'), grid on
    %xlim([-5,5]), ylim([5e-3,1]), xlabel('f'), %ylabel('Amplitude'), 
    %title(['DoF #' num2str(ii) ' in the Frequency Domain']), set(gca,'FontSize',13)
    drawnow
end

figure
for ii=1:4
    subplot(4,1,ii)
    semilogy(f,fftshift(squeeze(abs(Qhat(checkvec(end-1),:,ii*2))))/Nt,f,fftshift(squeeze(abs(Qhat_rom(checkvec(end-1),:,ii*2))))/Nt,'-.')
    legend('FOM','ROM')
    xlim([-20,20]), xlabel('f'), ylabel('Amplitude'), grid on
    title(['DoF #' num2str(checkvec(end-1)) ' at block ',num2str(ii*2),' in the Frequency Domain'])
    drawnow
end

figure
for ii=0:4
    subplot(5,1,ii+1)
    fi = 3*ii+2;
    plot(1:Nb,squeeze(abs(Aspod(1,:,fi))),1:Nb,squeeze(abs(A_rom(1,:,fi))),'-.')
    legend('Data','SPOD-ROM') , ylabel('Amp.'), grid on
    title(['A_{f' num2str(fi) '}']), 
    set(gca,'FontSize',9,'TickDir','out'), drawnow
end
    
%Perform inverse-DFT operation on every reconstructed block
Q_rom = zeros(Nx,Nf,Nb);
for n = 1:Nb
    Q_rom(:,:,n) = Nf/winWeight*ifft(squeeze(Qhat_rom(:,:,n)),Nf,2)./window; 
end
disp('Inverse DFT performed. SPOD_ROM obtained!')
Q_rom_real = real(Q_rom);

%{
imag_percentage = zeros(Nb,1);
for ii = 1:Nb
    imag_percentage(ii) = norm(imag(Q_rom(:,:,ii)))/norm(Q_rom_real(:,:,ii));
end
disp(['max imaginary percentage measured in 2-norm is ' num2str(max(imag_percentage))])
%}

figure;
ctr = 1;
for n = 1:5
    for tt = 1:(Nf-No)
        subplot(1,2,1)
        pcolor(xq,yq,reshape(Q(:,tt,n),[sqrt(Nx),sqrt(Nx)])); caxis([-3 3]), shading interp, axis equal tight,
        title({'FOM in block ' num2str(n), ' timestep tt=' num2str(tt), '#DoF=9604'})
        set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dx L_x-dx],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01])%,'XAxisLocation','top')
        drawnow
        
        subplot(1,2,2)
        pcolor(xq,yq,reshape(Q_rom_real(:,tt,n),[sqrt(Nx),sqrt(Nx)])); caxis([-3 3]), shading interp, axis equal tight, 
        title({'SPOD-ROM in block ' num2str(n), ' timestep tt=' num2str(tt), '#DoF=16 (4 SPOD modes at each f_k)'})
        set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dx L_x-dx],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01])%,'XAxisLocation','top')
        drawnow
        
        Film(ctr) = getframe(gcf);
        ctr = ctr+1;
    end
end

%{
figure;
ctr=1; tt_off = size(Xnoise,2)-size(X,2);
for tt = 205
    n = 1;
    subplot(1,3,1)
    pcolor(xq,yq,reshape(Q(:,tt,n),[sqrt(Nx),sqrt(Nx)])); caxis([-3 3]), shading interp, axis equal tight, colorbar
    %title(['HFM in block ' num2str(n) ', time = ' num2str(tt)])
    title('HFM')
    set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dx L_x-dx],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01])%,'XAxisLocation','top')
    drawnow
    %ctr = ctr + 1;

    subplot(1,3,2)
    pcolor(xq,yq,reshape(Q_rom_real(:,tt,n),[sqrt(Nx),sqrt(Nx)])); caxis([-3 3]), shading interp, axis equal tight, colorbar
    %title(['SPOD-ROM in block ' num2str(n), ', time = ' num2str(tt)])
    title(['SPOD-ROM b_k = ' num2str(bk)])
    set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dx L_x-dx],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01])%,'XAxisLocation','top')
    drawnow
    %ctr = ctr + 1;
    
    subplot(1,3,3)
    pcolor(xq,yq,reshape(X_rom(:,tt_off + tt),[sqrt(Nx),sqrt(Nx)])); caxis([-3 3]), shading interp, axis equal tight, colorbar
    %title(['POD-ROM in block ' num2str(n), ', time = ' num2str(tt)])
    title(['POD-ROM k = ' num2str(k)])
    set(gca,'Ydir','reverse','XTick',[dx L_x-dx],'XTickLabel',[0 1],'YTick',[dx L_x-dx],'YTickLabel',[0 1],'TickDir','out','TickLength',[0.02 0.01])%,'XAxisLocation','top')
    drawnow
end
%}

%% Compute errors
rel_errvec_fro = zeros(Nb,1);
rel_errvec_inf = zeros(Nb,1);
rel_errvec_2 = zeros(Nb,1);
for n = 1:Nb
    rel_errvec_fro(n,1) = norm(Q(:,:,n) - Q_rom_real(:,:,n),'fro')/norm(Q(:,:,n),'fro');
    rel_errvec_inf(n,1) = norm(Q(:,:,n) - Q_rom_real(:,:,n),'inf')/norm(Q(:,:,n),'inf');
    rel_errvec_2(n,1) = norm(Q(:,:,n) - Q_rom_real(:,:,n),2)/norm(Q(:,:,n),2);
end
mean(rel_errvec_fro)
mean(rel_errvec_inf)
mean(rel_errvec_2)

%{
rel_errvec_2 = zeros(1,Nb);
rel_err_2 = zeros(1,Nf);
for n = 1:Nb
    for tt = 1:Nf
        rel_err_2(tt) = norm(Q(:,tt,n) - Q_rom_real(:,tt,n),2)/norm(Q(:,tt,n),2)
    end 
    rel_errvec_2(n) = mean(rel_err_2);
end
mean(rel_errvec_2)
%}

%% Compute errors
rel_errvec_fro = zeros(Nb,Nf);
rel_errvec_inf = zeros(Nb,Nf);
rel_errvec_2 = zeros(Nb,Nf);
for n = 1:Nb
    for tt = 1:Nf
        rel_errvec_fro(n,tt) = norm(Q(:,tt,n) - Q_rom_real(:,tt,n),'fro')/norm(Q(:,tt,n),'fro');
        rel_errvec_inf(n,tt) = norm(Q(:,tt,n) - Q_rom_real(:,tt,n),'inf')/norm(Q(:,tt,n),'inf');
        rel_errvec_2(n,tt) = norm(Q(:,tt,n) - Q_rom_real(:,tt,n),2)/norm(Q(:,tt,n),2);
    end
end
figure
plot(rel_errvec_fro'), xlabel('timestep'), ylabel('error')
title('Relative frobenius-error over the realization time window')
figure
plot(rel_errvec_inf'), xlabel('timestep'), ylabel('error')
title('Relative infinity-error over the realization time window')
figure
plot(rel_errvec_2'), xlabel('timestep'), ylabel('error')
title('Relative 2-error over the realization time window')

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
