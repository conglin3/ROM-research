function a_dot = ode2DBurgersROM(a,Qrom,Lrom,Crom,Q3D,k)
%a_dot = Uk'*((Uk11*a).*(CXUk*a) + (Uk22*a).*(CYUk*a)) + Dstar*a;
a_dot = Qrom * kron(a,a) + Lrom*a + Crom; 
%{
%TENSOR FORM ROM
a_dot = zeros(k,1);
ajak = permute(a*a',[3,2,1]);
for ii = 1:k
   a_dot(ii) = sum(sum(Q3D(ii,:,:).*ajak)) + Lrom(ii,:)*a + Crom(ii); 
end
%}