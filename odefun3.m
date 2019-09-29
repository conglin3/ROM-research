function dydt = odefun3(y,A0star,Uk,Cx,u0,Qrom,Lrom,Crom)
   % the kxk matrix to solve for the time dependent weighing factors a_1(t) to a_k(t)
   % dydt = Uk'*(A0star + diag(Uk*y)*Cx)*(Uk*y + u0);
   dydt = Qrom * kron(y,y) + Lrom*y + Crom;     % precomputed form
end