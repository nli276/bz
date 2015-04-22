%% Tr,Yr
%solve ODEs (modified for light pulses)

% clear all
% close all

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
rtaumx=zeros(length(locs)-1,1);
for i=1:length(locs)-1
    rtaumx(i)=locs(i+1)-locs(i);
end
tau0=round(mean(rtaumx));

figure
plot(rtaumx)

icr=Yr(locs(1),:);
[Tr,Yr] = ode15s(@(t,y) turing(t,y,n,par),[0:1:runtime],icr,options);
figure
plot(Tr,Yr)

%% phi1dot=H(phi2-phi1) (bottom up), mathematically more rigorous as in the book, but essentially the same
% theta may be larger than 2pi
%calculate phidot numerically, need to compare this to eqn(10.15)
%Izhikevich using H calculated from run_turing.m
%dphi=dtheta, theta0=(1/tau0)*T.
n=2;%# of droplets -including bdy drops
% par=[0.05 0 0 0 2.5 1 0];%Px=0.05~0.1
par=[0 0 0 0 600 0 0];
%initial conditions for each droplet.
ic=zeros(1,7*n);

% ic(1)=X0; ic(2)=Y0; ic(7)=C0;%oil start from zero anyway
options = odeset('RelTol',1e-8,'AbsTol',1e-9);
%dudt not as good as 1e-12/1e-13, but Z does
ps=20+1;%step points
step=(locs(2)-locs(1))/(ps-1);%phase step size
tn=locs(2)-locs(1)+1;%period + 1, kinda like ps

theta1=zeros(runtime+1,ps,ps);%drop1
% theta01=zeros(runtime+1,ps);
phi1=zeros(runtime+1,ps,ps);
theta2=zeros(runtime+1,ps,ps);%drop2
% theta02=zeros(runtime+1,ps,ps);
phi2=zeros(runtime+1,ps,ps);
dphi=zeros(runtime+1,ps,ps);
% dtheta=zeros(runtime+1,ps,ps);
% dtheta0=zeros(runtime+1,ps,ps);
phasenorm1=zeros(tn,1);
phasenorm2=zeros(tn,1);
psi=[0:1/(ps-1):1];% initial difference, 1 is 2pi

%ref phase
theta0=zeros(runtime+1,1);%ref
for i=1:runtime+1
    theta0(i)=(i-1)/(tn-1);
end
figure
plot(theta0)
%%%%

for k=8:8%(ps-1)/2+1:(ps-1)/2+1%1:ps
    clear Y
    clear T
    clear I
    
    ostep=1+round((k-1)*step);%1~tn
    ic(1:7)=Yr(ostep,1:7);%drop 1
    %     ic(1:7)=icr(1:7);%icr is (second peak loc-step) ic of Yr(ic0)
    for i=1:ps
        nstep=round(ostep-(i-1)*step);%ostep~ostep-(tn-1)
        if nstep<0
            nstep=nstep+tn-1;
        end
        ic(8:14)=Yr(nstep,1:7);%drop 2
        [T,Y] = ode15s(@(t,y) turing(t,y,n,par),[0:1:runtime],ic,options);
        for j=1:runtime+1
            for l=1:tn
                phasenorm1(l)=norm(Y(j,1:7)-Yr(l,1:7));
                phasenorm2(l)=norm(Y(j,8:14)-Yr(l,1:7));
            end
            
            [C1,I1]=min(phasenorm1);
            [C2,I2]=min(phasenorm2);
            
            % deltheta(i,k)=(k-1)/(ps-1)-(I-1)/(tn-1);%H function
            theta1(j,i,k)=(I1-1)/(tn-1);
            if j>1
                while theta1(j,i,k)<theta1(j-1,i,k)-0.1
                    theta1(j,i,k)=theta1(j,i,k)+1;%theta keep increasing
                end
            end
            
            theta2(j,i,k)=(I2-1)/(tn-1);
            if j>1
                while theta2(j,i,k)<theta2(j-1,i,k)-0.1
                    theta2(j,i,k)=theta2(j,i,k)+1;
                end
            end
            
            phi1(j,i,k)=theta1(j,i,k)-theta0(j);
            phi2(j,i,k)=theta2(j,i,k)-theta0(j);
            
            dphi(j,i,k)=phi1(j,i,k)-phi2(j,i,k);
            %dphi=dtheta
        end
        
        if min(dphi(:,i,k))<0
            dphi(:,i,k)=dphi(:,i,k)+1;
        end
        
        if max(dphi(:,i,k))>1
            dphi(:,i,k)=dphi(:,i,k)-1;
        end
        figure(6)
        plot(Tr,Yr(:,s)+3e-3,'k')
        hold all
        plot(T,Y(:,s)+3e-3,'r',T,Y(:,s+7),'b')
        title(['\theta_1 #',num2str(k),'; \phi #',num2str(i)]);
        hold off
        
    end
end
%%
close all
% k=11;
i=8;

figure
plot(Tr,theta1(:,i,k),Tr,theta2(:,i,k),Tr,theta0)%Tr,theta01(:,k),'.',Tr,theta02(:,i,k),'.')
xlabel('time (s)')
ylabel('\theta_1,\theta_2, and reference in dots')

figure
% plot(Tr,phi1(:,i,k),Tr,phi2(:,i,k),Tr,dphi(:,:,k))
plot(Tr,dphi(:,:,k))
xlabel('time (s)')
% ylabel('\phi_1, \phi_2, \Delta\phi')
ylabel('\Delta\phi')

%% check eqn (10.15) Izhikevich (top down)
close all

Hchi=zeros(runtime+1,ps);%H1 in 10.3.1
Hchii=zeros(runtime+1,ps);%H2
psi1=zeros(runtime+1,ps);
psi2=zeros(runtime+1,ps);
psi1(1,:)=phi1(1,:,k);%

for i=3:4:19%1:ps%20:20%
    psi2(1,i)=psi1(1,i)-(i-1)/(ps-1);
    for j=1:runtime+1
        Hchi(j,i)=interp1(chi,H,dphi(j,i,k));
        Hchii(j,i)=interp1(chi,H,1-dphi(j,i,k));
        %0<dphi<1
        %chi=dtheta(t=0)=dphi(t=0)
        if j<runtime+1
            psi1(j+1,i)=psi1(j,i)+(1/tau0)*Hchi(j,i);
            psi2(j+1,i)=psi2(j,i)+(1/tau0)*Hchii(j,i);
            %epsilon=1, (1/tau0)*H is the actual H in (10.16)
            %integrate over time
        end
    end
    figure(1)
    hold all
    plot(chi,H,'k')
    plot(dphi(:,i,k),Hchi(:,i),dphi(:,i,k),Hchii(:,i),'linewidth',3)
    hold off
    
    figure(2)%compare psi with phi
%         hold all
    plot(Tr,phi1(:,i,k),Tr,psi1(:,i))
    xlabel('time (s)');
    ylabel('\phi_1 vs. \psi_1');
    title(['Initial phase difference=' num2str((i-1)/(ps-1))]);
%         hold off
    figure(3)
    hold all
    plot(Tr,dphi(:,i,k),Tr,(psi1(:,i)-psi2(:,i)))
    xlabel('time (s)');
    ylabel('\Delta\phi vs. \Delta\psi');
    hold off    
%     pause
end


