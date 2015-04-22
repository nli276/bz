%% Tr,Yr
%solve ODEs

clear all
close all

runtime=2000;%500;%
n=3;%# of droplets on a ring, n>=2

%%%%%%%%%%%%%%%chemistry%%%%%%%%%
h0=0.16;%[H+] Molar
A=0.3;%[BrO3-]
m=0.4;%[malonic acid]
%don't forget to change k9p if using low MA
%oscillation period increase as lowering m, so more runtime needed
%%%%%%%
cvh=0/100;%Randomness (the coefficient of variation) for acidity [H+], zero means monodispersed
h=abs(normrnd(h0,h0*0,1,n));
%%%%%%%%%%%%%%%geometry%%%%%%%%%%%%
lw = 0.01; %length of bz (water) drop(cm), 0.02 is 200um
cva=0/100; %The coefficient of variation (CV) equals the standard deviation divided by the mean (expressed as a percent).
lo = lw; %length of oil drop(cm)
cvb=0/100; %randomness in oil gap size

bz=abs(normrnd(lw,lw*cva,1,n));
% bz=sqrt(lw.*[0.0377    0.0035    0.0403]);%pi-s
oil=abs(normrnd(lo,lo*cvb,1,n));
% oil=sqrt(lw.*[0.0377    0.0035    0.0403]);

%%%%%%%%%%%%diffusion%%%%%%%%%%%%%%%
par=[0 0 0 0];%zero for ref partition coefficient for the 4 species, i.e. drops are not coupled for ref
dc=zeros(4,n);
D = 1e-5; %diffusion coefficient cm^2/s, I assume it's the same for all species
for i=1:4
    for j=1:n
        dc(i,j)=D*par(i)/(bz(j)*oil(j));
    end
end

ic=zeros(1,4*n);
%basic initial conditions
X0=1e-6;
Y0=5e-5;
Z0=0;
U0=0;

for i=1:n
    ic(1+(i-1)*4)=X0;
    ic(2+(i-1)*4)=Y0;
end
options = odeset('RelTol',1e-8,'AbsTol',1e-9);

% reference Tr, Yr
[Tr,Yr] = ode15s(@(t,y) fourvar(t,y,n,dc,h,A,m),[0:1:runtime],ic,options);
figure
plot(Tr,Yr)

% mph=1e-3;
mph=[1e-4 5e-4 1e-3 2e-5];%highMA(>=0.1M)
% mph=[1e-4 4e-5 1e-3 3e-5];%lowMA(<0.1M)
mpd=100;%not time, number of points in Y

% ref tau0
s=3;%Fe3+
clear pks
clear locs
[pks,locs]=findpeaks(Yr(:,s),'minpeakheight',mph(s),'minpeakdistance',mpd);
rtaumx=zeros(length(locs)-1,1);
for i=1:length(locs)-1
    rtaumx(i)=locs(i+1)-locs(i);
end
% tau0=round(mean(rtaumx));

figure
plot(rtaumx)
mtau=round(mean(rtaumx));

icr=Yr(locs(1),:);%star from zero phase

[Tr,Yr] = ode15s(@(t,y) fourvar(t,y,n,dc,h,A,m),[0:1:mtau],icr,options);
figure
plot(Tr,Yr)

%IC for every ten degree out of 360 (2pi)
ic360=zeros(1,4,36);
for j=1:36
    ic360(1,:,j)=Yr((j-1)*round(mtau/36)+1,1:4);
end
%% T,Y
%initial conditions (initial phase) from Tr,Yr
ips=randi(36,1,n);%random initial phase for each drop (10 degree resolution)
% ips=[27 12 2];
% ips=[3    32    23];

h=abs(normrnd(h0,h0*cvh,1,n));
% h=[0.1744 0.1645 0.1110];

for i=1:n
    ip(1+4*(i-1):4+4*(i-1))=ic360(1,1:4,ips(i));
end

par=[0.1 0 0 2.5];%partition coefficient ([oil]/[water] at interface)
%0.1 (estimated) for HBrO2 (weak acid), and 2.5 (measured) for Br2 (nonpolar)
dc=zeros(4,n);
D = 1e-5; %diffusion coefficient cm^2/s, I assume it's the same for both Br2 and radical
for i=1:4
    for j=1:n
        dc(i,j)=D*par(i)/(bz(j)*oil(j));
    end
end

[T,Y] = ode15s(@(t,y) fourvar(t,y,n,dc,h,A,m),[0:1:runtime],ip,options);

%%
Ys=zeros(size(Y));%shifted for visualization
%from bottom to top:Fe3+ (s=3) for drop 1, 2, 3, ...
figure
hold all
for i=1:n
    Ys(:,s+(i-1)*4)=Y(:,s+(i-1)*4)+(i-1)*3e-3;
    plot(T,Ys(:,s+(i-1)*4))
end
hold off
%% show sizes
bz
oil
h
% A
% m
ips