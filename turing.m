%created by Ning Li, based on eqn 6.1 in Turing's 1952 paper for diffusion
% and 7 variable FKN model for reaction
function dy = turing(t,y,n,par)
%simplified FKN model with constant BrO3-, MA, BrMA, and H+
%no flux boundaries

h = 0.08; %[H+](M)
A = 0.3; %[BrO3-]
m = 0.4; %[MA]
k1 = 2e6; %(1/Molsec)=(1/Ms)
k2 = 2; %(1/s)
k3 = 3000; %(1/Ms)
k4 = 42; %(1/s)
k5 = 5e9; %(1/Ms)
k6 = 10; %(1/s)
k7 = 29; %(1/s)
k8 = 9.3; %(1/s)
k9 = 0.1;%0.07;% for m>=0.1; 0.07 for m<0.1
b = k9*m;
k10 = 0.05; %(1/s)
kr = 2e8; %(1/Ms)
kred = 5e6; %(1/Ms)
bc=0.05; %M

%no light perturbation
ki=zeros(1,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 0.011; %length of water drop(cm)
lo = 0.008; %length of oil drop(cm)
D = 1e-5; %diffusion coefficient cm^2/s
dc = D*par./(lw*lo);%diffusion constant D*P(i)/ab

%y(1) = X, [HBrO2]
%y(2) = Y, [Br-]
%y(3) = Z, [Fe(III)] oxidized catalyst
%y(4) = p, [HOBr]
%y(5) = U, [Br2] in water
%y(6) = W, [BrO2*] in water
%y(7) = C, [Fe(II)] reduced catalyst


dy = zeros(7*n,1);
for i=0:n-1
    if i==0 %first drop (no flux)
        dy(1+i*7) = -k1*y(1+i*7)*y(2+i*7)*h+k2*y(2+i*7)*A*h^2-2*k3*y(1+i*7)^2-k4*y(1+i*7)*A*h+kr*y(6+i*7)^2+kred*y(6+i*7)*y(7+i*7) + dc(1)*(y(1+(i+1)*7)-y(1+i*7));
        dy(2+i*7) = -k1*y(1+i*7)*y(2+i*7)*h-k2*y(2+i*7)*A*h^2-k5*y(2+i*7)*y(4+i*7)*h+k6*y(5+i*7)+k7*y(5+i*7)*m+k9*y(3+i*7)*m+ki(i+1)*y(7+i*7)*b/(bc+b) + dc(2)*(y(2+(i+1)*7)-y(2+i*7));
        dy(3+i*7) = kred*y(6+i*7)*y(7+i*7)-k9*y(3+i*7)*m-k10*y(3+i*7)*m+ki(i+1)*y(7+i*7)*b/(bc+b) + dc(3)*(y(3+(i+1)*7)-y(3+i*7));
        dy(4+i*7) = 2*k1*y(1+i*7)*y(2+i*7)*h+k2*y(2+i*7)*A*h^2+k3*y(1+i*7)^2-k5*y(2+i*7)*y(4+i*7)*h+k6*y(5+i*7)-k8*y(4+i*7)*m + dc(4)*(y(4+(i+1)*7)-y(4+i*7));
        dy(5+i*7) = k5*y(2+i*7)*y(4+i*7)*h-k6*y(5+i*7)-k7*y(5+i*7)*m + dc(5)*(y(5+(i+1)*7)-y(5+i*7));
        dy(6+i*7) = 2*k4*y(1+i*7)*A*h-2*kr*y(6+i*7)^2-kred*y(6+i*7)*y(7+i*7) + dc(6)*(y(6+(i+1)*7)-y(6+i*7));
        dy(7+i*7) = -kred*y(6+i*7)*y(7+i*7)+k9*y(3+i*7)*m+k10*y(3+i*7)*m-ki(i+1)*y(7+i*7)*b/(bc+b) + dc(7)*(y(7+(i+1)*7)-y(7+i*7));
        
    elseif i==n-1 %last drop (no flux)
        dy(1+i*7) = -k1*y(1+i*7)*y(2+i*7)*h+k2*y(2+i*7)*A*h^2-2*k3*y(1+i*7)^2-k4*y(1+i*7)*A*h+kr*y(6+i*7)^2+kred*y(6+i*7)*y(7+i*7) + dc(1)*(y(1+(i-1)*7)-y(1+i*7));
        dy(2+i*7) = -k1*y(1+i*7)*y(2+i*7)*h-k2*y(2+i*7)*A*h^2-k5*y(2+i*7)*y(4+i*7)*h+k6*y(5+i*7)+k7*y(5+i*7)*m+k9*y(3+i*7)*m+ki(i+1)*y(7+i*7)*b/(bc+b) + dc(2)*(y(2+(i-1)*7)-y(2+i*7));
        dy(3+i*7) = kred*y(6+i*7)*y(7+i*7)-k9*y(3+i*7)*m-k10*y(3+i*7)*m+ki(i+1)*y(7+i*7)*b/(bc+b) + dc(3)*(y(3+(i-1)*7)-y(3+i*7));
        dy(4+i*7) = 2*k1*y(1+i*7)*y(2+i*7)*h+k2*y(2+i*7)*A*h^2+k3*y(1+i*7)^2-k5*y(2+i*7)*y(4+i*7)*h+k6*y(5+i*7)-k8*y(4+i*7)*m + dc(4)*(y(4+(i-1)*7)-y(4+i*7));
        dy(5+i*7) = k5*y(2+i*7)*y(4+i*7)*h-k6*y(5+i*7)-k7*y(5+i*7)*m + dc(5)*(y(5+(i-1)*7)-y(5+i*7));
        dy(6+i*7) = 2*k4*y(1+i*7)*A*h-2*kr*y(6+i*7)^2-kred*y(6+i*7)*y(7+i*7) + dc(6)*(y(6+(i-1)*7)-y(6+i*7));
        dy(7+i*7) = -kred*y(6+i*7)*y(7+i*7)+k9*y(3+i*7)*m+k10*y(3+i*7)*m-ki(i+1)*y(7+i*7)*b/(bc+b) + dc(7)*(y(7+(i-1)*7)-y(7+i*7));
    
    else %drops in between
        dy(1+i*7) = -k1*y(1+i*7)*y(2+i*7)*h+k2*y(2+i*7)*A*h^2-2*k3*y(1+i*7)^2-k4*y(1+i*7)*A*h+kr*y(6+i*7)^2+kred*y(6+i*7)*y(7+i*7) + dc(1)*(y(1+(i+1)*7)+y(1+(i-1)*7)-2*y(1+i*7));
        dy(2+i*7) = -k1*y(1+i*7)*y(2+i*7)*h-k2*y(2+i*7)*A*h^2-k5*y(2+i*7)*y(4+i*7)*h+k6*y(5+i*7)+k7*y(5+i*7)*m+k9*y(3+i*7)*m+ki(i+1)*y(7+i*7)*b/(bc+b) + dc(2)*(y(2+(i+1)*7)+y(2+(i-1)*7)-2*y(2+i*7));
        dy(3+i*7) = kred*y(6+i*7)*y(7+i*7)-k9*y(3+i*7)*m-k10*y(3+i*7)*m+ki(i+1)*y(7+i*7)*b/(bc+b) + dc(3)*(y(3+(i+1)*7)+y(3+(i-1)*7)-2*y(3+i*7));
        dy(4+i*7) = 2*k1*y(1+i*7)*y(2+i*7)*h+k2*y(2+i*7)*A*h^2+k3*y(1+i*7)^2-k5*y(2+i*7)*y(4+i*7)*h+k6*y(5+i*7)-k8*y(4+i*7)*m + dc(4)*(y(4+(i+1)*7)+y(4+(i-1)*7)-2*y(4+i*7));
        dy(5+i*7) = k5*y(2+i*7)*y(4+i*7)*h-k6*y(5+i*7)-k7*y(5+i*7)*m + dc(5)*(y(5+(i+1)*7)+y(5+(i-1)*7)-2*y(5+i*7));
        dy(6+i*7) = 2*k4*y(1+i*7)*A*h-2*kr*y(6+i*7)^2-kred*y(6+i*7)*y(7+i*7) + dc(6)*(y(6+(i+1)*7)+y(6+(i-1)*7)-2*y(6+i*7));
        dy(7+i*7) = -kred*y(6+i*7)*y(7+i*7)+k9*y(3+i*7)*m+k10*y(3+i*7)*m-ki(i+1)*y(7+i*7)*b/(bc+b) + dc(7)*(y(7+(i+1)*7)+y(7+(i-1)*7)-2*y(7+i*7));
        
    end;
end;
