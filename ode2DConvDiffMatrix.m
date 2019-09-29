function u_dot = ode2DConvDiffMatrix(t,y,Asys,bsys,Nx,Ny,r1,r2,f1,f2,dt,noise)
   % the FOM
   tt = round(t/dt)+1;
   u_dot = Asys*y + bsys + r1*cos(2*pi*f1*t) + r2*cos(2*pi*f2*t) + noise(:,tt);
end