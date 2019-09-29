function a_dot = ode2DConvDiffPODROM(t,y,A_rom,b_rom,c_rom,Nx,Ny,r1_rom,r2_rom,f1,f2,dt,noise_rom)
   % the FOM
   tt = round(t/dt)+1;
   a_dot = A_rom*y + b_rom + c_rom + r1_rom*cos(2*pi*f1*t) + r2_rom*cos(2*pi*f2*t) + noise_rom(:,tt);
end