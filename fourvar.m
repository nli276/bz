%created by Ning Li on Mar 28, 2012, based on eqn 6.1 in Turing's 1952 paper
% and eqn (29)-(32) in ref 2 (Vanag & Epstein, 2009)
%Periodic boundary condition and No flux boundary condition
%no light perturbation
%drop size and chemical heterogeneity is available
function dy = fourvar(t,y,n,dc,h,A,m)

% h = 0.16; %[H+](Mole)
% A = 0.3; %[BrO3-]
% m = 0.4; %[MA]

k1 = 2e6*h; %(1/Molsec)=(1/Ms)
k2 = 2*A*h.^2; %(1/s)
k3 = 3000; %(1/Ms)
k4 = 42*A*h; %(1/s)
k7 = 29*m; %(1/s)
k9p=0.1;%0.1 for m>+0.1;0.07 for m<0.1;% (1/s),k9p=0~0.2
k9 = k9p*m;
k10 = 0.05*m; %(1/s)
kr = 2e8; %(1/Ms)
kred = 5e6; %(1/Ms)

c0=3e-3;
cmin=sqrt(2*kr*(k9+k10)*c0/kred^2);
%%%%%%%%%%%%%%%%% four-variable ODEs %%%%%%%%%%%%%%%%%%%
%y(1) = X, [HBrO2]
%y(2) = Y, [Br-]
%y(3) = Z, [Fe(III)] oxidized catalyst
%y(4) = U, [Br2]


dy = zeros(4*n,1);
%%%% for Periodic boundary condition (drops on a ring) %%%%%%%%%
for i=0:n-1
    if i==0 %first drop (pbc)
        dy(1+i*4) = -k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)-2*k3*y(1+i*4)^2+k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin) + dc(1,i+1)*(y(1+(i+1)*4)+y(1+(n-1)*4)-2*y(1+i*4));
        dy(2+i*4) = -3*k1(i+1)*y(1+i*4)*y(2+i*4)-2*k2(i+1)*y(2+i*4)-k3*y(1+i*4)^2+k7*y(4+i*4)+k9*y(3+i*4)+ dc(2,i+1)*(y(2+(i+1)*4)+y(2+(n-1)*4)-2*y(2+i*4));
        dy(3+i*4) = 2*k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)-k9*y(3+i*4)-k10*y(3+i*4)+ dc(3,i+1)*(y(3+(i+1)*4)+y(3+(n-1)*4)-2*y(3+i*4));%-0.001*y(3+i*4);%ferroin decay?
        dy(4+i*4) = 2*k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)+k3*y(1+i*4)^2-k7*y(4+i*4)+ dc(4,i+1)*(y(4+(i+1)*4)+y(4+(n-1)*4)-2*y(4+i*4));

    elseif i==n-1 %last drop (pbc)
        dy(1+i*4) = -k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)-2*k3*y(1+i*4)^2+k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)+ dc(1,i+1)*(y(1+0*4)+y(1+(i-1)*4)-2*y(1+i*4));
        dy(2+i*4) = -3*k1(i+1)*y(1+i*4)*y(2+i*4)-2*k2(i+1)*y(2+i*4)-k3*y(1+i*4)^2+k7*y(4+i*4)+k9*y(3+i*4)+ dc(2,i+1)*(y(2+0*4)+y(2+(i-1)*4)-2*y(2+i*4));
        dy(3+i*4) = 2*k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)-k9*y(3+i*4)-k10*y(3+i*4)+ dc(3,i+1)*(y(3+0*4)+y(3+(i-1)*4)-2*y(3+i*4));%-0.001*y(3+i*4);%ferroin decay?
        dy(4+i*4) = 2*k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)+k3*y(1+i*4)^2-k7*y(4+i*4)+ dc(4,i+1)*(y(4+0*4)+y(4+(i-1)*4)-2*y(4+i*4));
       
    else %drops in between
        dy(1+i*4) = -k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)-2*k3*y(1+i*4)^2+k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin) + dc(1,i+1)*(y(1+(i+1)*4)+y(1+(i-1)*4)-2*y(1+i*4));
        dy(2+i*4) = -3*k1(i+1)*y(1+i*4)*y(2+i*4)-2*k2(i+1)*y(2+i*4)-k3*y(1+i*4)^2+k7*y(4+i*4)+k9*y(3+i*4)+ dc(2,i+1)*(y(2+(i+1)*4)+y(2+(i-1)*4)-2*y(2+i*4));
        dy(3+i*4) = 2*k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)-k9*y(3+i*4)-k10*y(3+i*4)+ dc(3,i+1)*(y(3+(i+1)*4)+y(3+(i-1)*4)-2*y(3+i*4));%-0.001*y(3+i*4);%ferroin decay?
        dy(4+i*4) = 2*k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)+k3*y(1+i*4)^2-k7*y(4+i*4)+ dc(4,i+1)*(y(4+(i+1)*4)+y(4+(i-1)*4)-2*y(4+i*4));
        
    end;
end;

%%%%%% No flux boundary condition %%%%%%%%%%%%%%%%%%%
% % 
% for i=0:n-1
%     if i==0 %first drop (no flux)
%         dy(1+i*4) = -k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)-2*k3*y(1+i*4)^2+k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)+ dc(1)*(y(1+(i+1)*4)-y(1+i*4));
%         dy(2+i*4) = -3*k1(i+1)*y(1+i*4)*y(2+i*4)-2*k2(i+1)*y(2+i*4)-k3*y(1+i*4)^2+k7*y(4+i*4)+k9*y(3+i*4)+ dc(2)*(y(2+(i+1)*4)-y(2+i*4));
%         dy(3+i*4) = 2*k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)-k9*y(3+i*4)-k10*y(3+i*4) + dc(3)*(y(3+(i+1)*4)-y(3+i*4));
%         dy(4+i*4) = 2*k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)+k3*y(1+i*4)^2-k7*y(4+i*4) + dc(4)*(y(4+(i+1)*4)-y(4+i*4));
%         
%     elseif i==n-1 %last drop (no flux)
%         dy(1+i*4) = -k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)-2*k3*y(1+i*4)^2+k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)+ dc(1)*(y(1+(i-1)*4)-y(1+i*4));
%         dy(2+i*4) = -3*k1(i+1)*y(1+i*4)*y(2+i*4)-2*k2(i+1)*y(2+i*4)-k3*y(1+i*4)^2+k7*y(4+i*4)+k9*y(3+i*4)+ dc(2)*(y(2+(i-1)*4)-y(2+i*4));
%         dy(3+i*4) = 2*k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)-k9*y(3+i*4)-k10*y(3+i*4) + dc(3)*(y(3+(i-1)*4)-y(3+i*4));
%         dy(4+i*4) = 2*k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)+k3*y(1+i*4)^2-k7*y(4+i*4) + dc(4)*(y(4+(i-1)*4)-y(4+i*4));
%        
%     else %drops in between
%         dy(1+i*4) = -k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)-2*k3*y(1+i*4)^2+k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)+ dc(1)*(y(1+(i+1)*4)+y(1+(i-1)*4)-2*y(1+i*4));
%         dy(2+i*4) = -3*k1(i+1)*y(1+i*4)*y(2+i*4)-2*k2(i+1)*y(2+i*4)-k3*y(1+i*4)^2+k7*y(4+i*4)+k9*y(3+i*4)+ dc(2)*(y(2+(i+1)*4)+y(2+(i-1)*4)-2*y(2+i*4));
%         dy(3+i*4) = 2*k4(i+1)*y(1+i*4)*(c0-y(3+i*4))/(c0-y(3+i*4)+cmin)-k9*y(3+i*4)-k10*y(3+i*4) + dc(3)*(y(3+(i+1)*4)+y(3+(i-1)*4)-2*y(3+i*4));
%         dy(4+i*4) = 2*k1(i+1)*y(1+i*4)*y(2+i*4)+k2(i+1)*y(2+i*4)+k3*y(1+i*4)^2-k7*y(4+i*4) + dc(4)*(y(4+(i+1)*4)+y(4+(i-1)*4)-2*y(4+i*4));
%         
%     end;
% end;


