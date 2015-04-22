%created by Ning Li on 2-22, 2011.Modified for light pulses
%Based on Vladimir's paper (J. Chem. Phys. 131, 2009) for reaction with light
%and Marcin's model for diffusion--coupling.

%there is strength and time control of light perturbation
function dy = fknode(t,y,n,lp,taop,duty,tl)
%simplified FKN model with constant BrO3-, MA, BrMA, and H+
%seven variable version with light option
%symmetric boundary oil droplets
% n=4;%number of droplets (always two extra for boundaries)
h = 0.16; %[H+](Mole)
A = 0.3; %[BrO3-]
m = 0.2; %[MA]
b = 0.1*m; %[BrMA]
k1 = 2e6; %(1/Molsec)=(1/Ms)
k2 = 2; %(1/s)
k3 = 3000; %(1/Ms)
k4 = 42; %(1/s)
k5 = 5e9; %(1/Ms)
k6 = 10; %(1/s)
k7 = 29; %(1/s)
k8 = 9.3; %(1/s)
k9 = 0.1; %(1/s),k9''
k10 = 0.05; %(1/s)
kr = 2e8; %(1/Ms)
kred = 5e6; %(1/Ms)
bc=0.05; %M

%light

%%%%%light pulses%%%%%
%taop: period of light pulses

%%%------bdy drops-----always on
% lpb=1e-2;%5e-3;%1e-4;%
% lpb=0;
% ki(1)=lpb;
% ki(n)=lpb;
%%%------drops between bdy-----
%here tl is the time to start light pulse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=2:n-1 %drop 1 perturbed
for i=1:n%drop 1 perturbed-no ltbdy
%     if i==2;
    if i==1
%         ki(i)=lp*((1+sign(t-tl))*0.5*(square(2*pi*(t-tl)/taop,duty)+1)/2); %light pulses
        ki(i)=lp*(1-(1+sign(t-tl))*0.5*(square(2*pi*(t-tl)/taop,duty)+1)/2); %light gates
%         ki(i)=lp*(1-sign(t-tl))*0.5;%light switch
    else
        ki(i)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=2:2:n-1 %light every other one for odd n
%     ki(i)=lp*((1+sign(t-tl))*0.5*(square(2*pi*(t-tl)/taop,duty)+1)/2); %light pulses
%     %ki(i)=lp*(1-(1+sign(t-tl))*0.5*(square(2*pi*(t-tl)/taop,duty)+1)/2); %light gates
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 0.02; %length water drop(cm)a~200um
lo = 0.01; %length oil drop(cm)b~50um
D = 1e-5; %diffusion coefficient cm^2/s, I assume it's the same for both Br2 and radical
pB = 2.5;%2.5;% partition coefficient for [Br2], Co*/Cw* at equilibrium,
pR = 1;%1;% partition coefficient for [BrO2*], NOT sure about the value:(
sl = lw + lo;% sum of water drop length lw and oil drop length lo
cdw = 2*D/(lw*sl);%effective reaction rate due to diffusion, water to oil
cdo = 2*D/(lo*sl);%effective reaction rate due to diffusion, oil to water
%y(1) = X, [HBrO2]
%y(2) = Y, [Br-]
%y(3) = Z, [Fe(III)] oxidized catalyst
%y(4) = p, [HOBr]
%y(5) = U, [Br2] in water
%y(6) = W, [BrO2*] in water
%y(7) = C, [Fe(II)] reduced catalyst

%y(8) = UU=S, [Br2] in oil to the right
%y(9) = WW=R, [BrO2*] in oil to the right

dy = zeros(9*n+2,1);%The output must be a column vector
for i=0:n-1
    dy(1+i*9) = -k1*y(1+i*9)*y(2+i*9)*h+k2*y(2+i*9)*A*h^2-2*k3*y(1+i*9)^2-k4*y(1+i*9)*A*h+kr*y(6+i*9)^2+kred*y(6+i*9)*y(7+i*9);
    dy(2+i*9) = -k1*y(1+i*9)*y(2+i*9)*h-k2*y(2+i*9)*A*h^2-k5*y(2+i*9)*y(4+i*9)*h+k6*y(5+i*9)+k7*y(5+i*9)*m+k9*y(3+i*9)*m+ki(i+1)*y(7+i*9)*b/(bc+b);
    dy(3+i*9) = kred*y(6+i*9)*y(7+i*9)-k9*y(3+i*9)*m-k10*y(3+i*9)*m+ki(i+1)*y(7+i*9)*b/(bc+b);%-0.001*y(3+i*9);%ferroin decay?
    dy(4+i*9) = 2*k1*y(1+i*9)*y(2+i*9)*h+k2*y(2+i*9)*A*h^2+k3*y(1+i*9)^2-k5*y(2+i*9)*y(4+i*9)*h+k6*y(5+i*9)-k8*y(4+i*9)*m;
    dy(7+i*9) = -kred*y(6+i*9)*y(7+i*9)+k9*y(3+i*9)*m+k10*y(3+i*9)*m-ki(i+1)*y(7+i*9)*b/(bc+b);
    
    if i==0
        dy(5+i*9) = k5*y(2+i*9)*y(4+i*9)*h-k6*y(5+i*9)-k7*y(5+i*9)*m + cdw*(y(8+i*9)+y(9*n+1)-2*pB*y(5+i*9));
        dy(6+i*9) = 2*k4*y(1+i*9)*A*h-2*kr*y(6+i*9)^2-kred*y(6+i*9)*y(7+i*9) + cdw*(y(9+i*9)+y(9*n+2)-2*pR*y(6+i*9));
    else
        dy(5+i*9) = k5*y(2+i*9)*y(4+i*9)*h-k6*y(5+i*9)-k7*y(5+i*9)*m + cdw*(y(8+(i-1)*9)+y(8+i*9)-2*pB*y(5+i*9));
        dy(6+i*9) = 2*k4*y(1+i*9)*A*h-2*kr*y(6+i*9)^2-kred*y(6+i*9)*y(7+i*9) + cdw*(y(9+(i-1)*9)+y(9+i*9)-2*pR*y(6+i*9));
    end;
    
    if i==n-1
        dy(8+i*9) =  cdo*(pB*y(5+i*9) - y(8+i*9));
        dy(9+i*9) =  cdo*(pR*y(6+i*9) - y(9+i*9));
    else
        dy(8+i*9) =  cdo*(pB*(y(5+i*9)+y(5+(i+1)*9)) - 2*y(8+i*9));
        dy(9+i*9) =  cdo*(pR*(y(6+i*9)+y(6+(i+1)*9)) - 2*y(9+i*9));
    end;
    
    dy(9*n+1) = cdo*(pB*y(5) - y(9*n+1));%[Br2] in the oil drop before drop1 to make bdy symmetric
    dy(9*n+2) = cdo*(pR*y(6) - y(9*n+2));%[BrO2*] in the oil drop before drop1 to make bdy symmetric
    
    
end;
