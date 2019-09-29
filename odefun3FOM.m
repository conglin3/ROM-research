function dydt = odefun3FOM(y,Cx,Dx)
   dydt = Dx*y + y.*(Cx*y);        %non conservative form
   %dydt = Dx*y + 0.5*Cx*(y.*y);     %conservative form
end