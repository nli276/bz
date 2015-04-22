%% Tr,Yr
%solve ODEs

clear all
close all

n=2;%# of droplets -including bdy drops
par=zeros(1,7);%for ref
runtime=2000;%500;%

%basic initial conditions
X0=1e-6;
Y0=5e-5;
Z0=0;
P0=0;
U0=0;
W0=0;
C0=0.003;


%initial conditions for each droplet.
ic=zeros(1,7*n);
for i=1:n
ic(1+(i-1)*7)=X0; 
ic(2+(i-1)*7)=Y0; 
ic(7+(i-1)*7)=C0;
end
options = odeset('RelTol',1e-8,'AbsTol',1e-9);

% reference Tr, Yr
[Tr,Yr] = ode15s(@(t,y) turing(t,y,n,par),[0:1:runtime],ic,options);
figure
plot(Tr,Yr)

% mph=1e-3;
% mph=[1e-4 5e-4 1e-3 5e-7 2e-5 1e-7 1e-3];%highMA
mph=[1e-4 4e-5 2.5e-3 2e-4 3e-5 5e-7 1e-3];%lowMA
mpd=100;%not time, number of points in Y

% ref tau0
s=3;%Ferroin
clear pks
clear locs
[pks,locs]=findpeaks(Yr(:,s),'minpeakheight',mph(s),'minpeakdistance',mpd);
% IS=zeros(runtime+1);
% IS=0.5*(Yr(:,2)+Yr(:,3).*3^2+Yr(:,7).*2^2+0.16);
% (max(IS)-min(IS))/mean(IS)
% plot(Tr,IS)
% ylabel('ionic strength')
rtaumx=zeros(length(locs)-1,1);
for i=1:length(locs)-1
    rtaumx(i)=locs(i+1)-locs(i);
end
tau0=round(mean(rtaumx))

figure
plot(rtaumx)

icr=Yr(locs(1),:);
[Tr,Yr] = ode15s(@(t,y) turing(t,y,n,par),[0:1:runtime],icr,options);
figure
plot(Tr,Yr)

Yrmax=max(Yr);
%% T,Y
n=2;%# of droplets -including bdy drops
% par=[0.1 0 0 0 2.5 0 0];
par=[0 0 0 0 0 0 0];%
%initial conditions for each droplet.
ic=zeros(1,7*n);

options = odeset('RelTol',1e-8,'AbsTol',1e-9);
ps=10+1;%step points
step=(locs(2)-locs(1))/(ps-1);%phase step size
tn=locs(2)-locs(1)+1;%period+1
phasenorm=zeros(tn,1);
phasenormwt=zeros(tn,1);
deltheta=zeros(ps,ps);
phi1=zeros(ps,ps);
phi1w=zeros(ps,ps);
phi=[0:1/(ps-1):1];% 1 is 2pi

for k=6:6%(ps-1)/2+1:(ps-1)/2+1%1:ps%1:ps%
    clear Y
    clear T
    clear I
    
    ostep=1+round((k-1)*step);%1~tn
    ic(1:7)=Yr(ostep,1:7);%drop 1
    for i=1:ps%11:11%3:3%
        nstep=round(ostep-(i-1)*step);%ostep~ostep-(tn-1)
        if nstep<0
            nstep=nstep+tn-1;
        end
        ic(8:14)=Yr(nstep,1:7);%drop 2
        [T,Y] = ode15s(@(t,y) turing(t,y,n,par),[0:1:5*tau0],ic,options);
        for j=1:tn
            phasenorm(j)=norm(Y(tn,1:7)-Yr(j,1:7));
            phasenormwt(j)=norm((Y(tn,1:7)-Yr(j,1:7))./(Yrmax(1:7)));
        end
        [C,I]=min(phasenorm);
%         deltheta(i,k)=(k-1)/(ps-1)-(I-1)/(tn-1);
        %phi1, actually it's ...+1-(...+1)
        phi1(i,k)=(I-1)/(tn-1)-(k-1)/(ps-1);
        %this is H*tau0
        
        %weighted
        [Cw,Iw]=min(phasenormwt);
        phi1w(i,k)=(Iw-1)/(tn-1)-(k-1)/(ps-1);
        
        figure(6)
        plot(Tr,Yr(:,s)+3e-3,'k')
        hold all
        plot(T,Y(:,s)+3e-3,'r',T,Y(:,s+7),'b')
        title(['\theta_1 #',num2str(k),'; \phi #',num2str(i)]);
        hold off
%         pause
    end
end
%%
% H=deltheta(:,k);
H=phi1(:,k);
Hw=phi1w(:,k);
chi=transpose(phi);
figure
plot(chi,H,chi,Hw)
title(['H function, [MA]=200 mM, par = ' num2str(par)]);
ylabel({'\Delta\theta'});
xlabel({'\theta'});
% figure
% plot(T,Y(:,s)+3e-3,'r',T,Y(:,s+7),'b')
% ylabel({'Ferroin (M)'});
% xlabel({'Time (s)'});