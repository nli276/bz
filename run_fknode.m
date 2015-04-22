%% solve ODEs (modified for light pulses)
clear
n=2;%# of droplets -including bdy drops
tn=1;%number of different ptb trials
% nt=1000;
runtime=1000;

%basic initial conditions
X0=1e-6;
Y0=5e-5;%increased Br- to 5e-5 get rid of the first "sync" peak
Z0=0;
P0=0;
U0=0;
W0=0;
C0=0.003;
UU0=0;
WW0=0;
%A0=0.3;
%M0=0.6;
%B0=0.2;
%H0=0.16;


%initial conditions for each droplet.
ic=zeros(1,9*n+2);
for i=0:n-1
    ic(1+i*9)=X0;
    ic(2+i*9)=Y0;
    ic(7+i*9)=C0;
end;

options = odeset('RelTol',1e-8,'AbsTol',1e-9);
%dudt not as good as 1e-12/1e-13, but Z does

%%%%----light parameters for fknode---loop one at a time for the easier coding--
% lp=2e-2;%5e-3;%1e-4;%%
lp=0;
step=200;%20;%
% taop1=10;
% taop=taop1:step:taop1+step*(tn-1);
taop=1e4;%1e3*2;%
% pud=10;%# seconds for light pulse duration
% duty=pud*(1./taop)*100;%percentage nubmer, 20 is 20% on
duty=18;%constant fraction of light

tl1=1000;
tl=tl1:step:tl1+step*(tn-1);
% tl=tl1;

pn=4;%peak number
mn=n;%# of droplets between silent boundary droplets
pks7=zeros(pn,mn*7,tn);%3 means z, we only care about that
locs7=zeros(pn,mn*7,tn);
tim7=zeros(pn,mn,tn);
% mph=1e-3;
% mph=[1e-4 5e-4 1e-3 1e-6 2e-5 5e-7 1e-3];%hightMA
mph=[1e-4 3.7e-5 1e-3 2e-4 2e-5 5e-7 1e-3];%lowMA
% mph=[1e-8 1e-8 1e-8 1e-8 1e-8 1e-8 1e-8];
mpd=50;%not time, number of points in Y

for i=1:tn
    clear Y
    clear T
    clear lt
    %         [T,Y] = ode15s(@(t,y) fknode(t,y,n,lp,taop(i),duty,tl),[0 runtime],ic,options);
    [T,Y] = ode15s(@(t,y) fknode(t,y,n,lp,taop,duty,tl(i)),[0 runtime],ic,options);%trick from Micha
    
    
    %%%%%%% show light %%%%%%%%%
    for k=1:length(T)
        %       lt(k)=lp*((1+sign(T(k)-tl))*0.5*(square(2*pi*(T(k)-tl)/taop(i),duty)+1)/2);
        %       lt(k)=lp*((1+sign(T(k)-tl(i)))*0.5*(square(2*pi*(T(k)-tl(i))/taop,duty)+1)/2);
                lt(k)=lp*(1-(1+sign(T(k)-tl(i)))*0.5*(square(2*pi*(T(k)-tl(i))/taop,duty)+1)/2);
%         lt(k)=lp*(1-sign(T(k)-tl(i)))*0.5;%light switch
    end
    %%%%%% find peaks for all 7 species (FKN->Kuramoto!) %%%%
    clear pks
    clear locs
    for d=1:mn%drop number, 1 is the first drop next to bdy
        if d==1%for one BZ pulse gated by light
            for s=1:7
%                 [pks,locs]=findpeaks(Y(:,s+9*d),'minpeakdistance',mpd);%'minpeakheight',mph(s),
                [pks,locs]=findpeaks(Y(:,s+9*(d-1)),'minpeakheight',mph(s),'minpeakdistance',mpd);%
                pks7(1,s,i)=pks(1);
                locs7(1,s,i)=locs(1);
                tim7(1,s,i)=T(locs(1));%7 species have different peak time
            end
            figure(d+mn)% show concentrations and light
            %plot(T,Y(:,1+9*d:7+9*d),T,lt)
            plot(T,Y(:,1+9*(d-1):7+9*(d-1)),T,lt)

            %title(['Tp=',num2str(taop(i)),'s']);
            title(['Start light at=',num2str(tl(i)),'s']);
        else
            for s=1:7
%                 [pks,locs]=findpeaks(Y(:,s+9*d),'minpeakheight',mph(s),'minpeakdistance',mpd);%
                [pks,locs]=findpeaks(Y(:,s+9*(d-1)),'minpeakheight',mph(s),'minpeakdistance',mpd);
                pks7(1:pn,s+7*(d-1),i)=pks(1:pn);
                locs7(1:pn,s+7*(d-1),i)=locs(1:pn);
                for j=1:pn
                    tim7(j,s+7*(d-1),i)=T(locs(j));%it's nonlinear, get the right time!
                end
            end
            figure(d+mn)% show concentrations and light
%             plot(T,Y(:,1+9*d:7+9*d))
            plot(T,Y(:,1+9*(d-1):7+9*(d-1)))
            %title(['Tp=',num2str(taop(i)),'s']);
            title(['Start light at=',num2str(tl(i)),'s']);
            
        end
        
        figure (d)
%         plot(T,Y(:,1+9*d:7+9*d))
        plot(T,Y(:,1+9*(d-1):7+9*(d-1)))
        hold all
        plot(tim7(:,1+7*(d-1):7*d,i),pks7(:,1+7*(d-1):7*d,i), 'marker', 'o');
        title(['Drop',num2str(d)]);
        hold off
        
    end
    
end

%% show it
ny=3;%select the chemical want to see
for j=ny:9:ny+9*(n-1)
    figure (j)
    plot(T,Y(:,j))
end

%% calculate tau0: manually produce a ref first
% pn=20;
% d=2;
rtaumx=zeros(pn-1,7);
for j=1:7
    for i=1:pn-1
        rtaumx(i,j)=tim7r(i+1,j+7,1)-tim7r(i,j+7,1);
    end
end
plot(rtaumx)
tau0=mean(rtaumx);%avg for each species.

%% calculate thetar (ref)
thetar=zeros(runtime+1,7*mn,tn);
time(:,1)=(0:runtime);%instead of taking from T, create smooth time
tim7ro=round(tim7r);
for i=1:tn
    for d=2:mn
        for j=1:pn-1
            for s=1:7
                ki=tim7ro(j,s+7*(d-1),i);
                kf=tim7ro(j+1,s+7*(d-1),i);
                for k=ki:kf-1
                    thetar(k+2,s+7*(d-1),i)=(k+1-ki)*2*pi/(kf-ki)+2*pi*(j-1);
                end
            end
        end
    end
    
end
%% plot thetar
figure(100)
plot(time,thetar(:,8:14,7))



%% calculate theta
theta=zeros(runtime+1,7*mn,tn);
time(:,1)=(0:runtime);%instead of taking from T, create smooth time
tim7o=round(tim7);
for i=1:tn
    for d=2:mn
        for j=1:pn-1
            for s=1:7
                ki=tim7o(j,s+7*(d-1),i);
                kf=tim7o(j+1,s+7*(d-1),i);
                for k=ki:kf-1
                    theta(k+2,s+7*(d-1),i)=(k+1-ki)*2*pi/(kf-ki)+2*pi*(j-1);
                end
            end
        end
    end
    
end
%% plot theta
figure(107)
hold all
% for i=1:tn
plot(time,theta(:,8:14,7))
% end
hold off
%% calculate instantaneous f(theta2-theta1)=dtheta/dt-w0; dt=1 as the stepsize in "time"
% this is just a number at the moment of drop1's peak that should fall on
% the next f curve. This is trivial.

theta21=zeros(tn,7);%drop1 has 1 pulse only
dtheta=zeros(tn,7);
fi=zeros(tn,7);%f instantaneous
for i=1:tn
    for d=2:mn
        for s=1:7
            k=tim7o(1,s,i);
%             theta21(i,s)=theta(k+1,s+7*(d-1),i);%this is x-axis
            theta21(i,s)=thetar(k+1,s+7*(d-1),i);
            dtheta(i,s)=theta(k+1,s+7*(d-1),i)-theta(k,s+7*(d-1),i);
            w0=2*pi./tau0;
            fi(i,s)=dtheta(i,s)-w0(s);
        end
    end
end
figure
plot(theta21,fi)

%% calculate f(theta21)=deltheta/tau0-w0;
% f is a function of phase, not time
% can use theta12, w0 from previous cell

% theta21=zeros(tn,7);
deltheta=zeros(tn,7);
f=zeros(tn,7);
j=1;
for i=1:tn
    for d=2:mn
        for s=1:7
            %           k=tim7o(1,s,i);
            %           theta21(i,s)=theta(k+1,s+7*(d-1),i);%this is x-axis
            %           w0=2*pi./tau0;
            
            while tim7ro(j,s+7*(d-1),i)<=tim7o(1,s,i)
                if j>pn-1, break, end
                j=j+1;
            end
            k0=tim7ro(j,s+7*(d-1),i);%from ref, so deltheta is absolute
            deltheta(i,s)=theta(k0+1,s+7*(d-1),i);
            f(i,s)=deltheta(i,s)./tau0(s)-(j-1)*w0(s);
        end
    end
end
figure
plot(theta21,f)

%% df

df=zeros(tn-1,7);
for i=1:tn-1
    for s=1:7
%         df(i,s)=(fi(i+1,s)-fi(i,s))/(theta21(i+1,s)-theta21(i,s));
        df(i,s)=(f(i+1,s)-f(i,s))/(theta21(i+1,s)-theta21(i,s));
    end
end
figure
plot(theta21(1:tn-1,:),df)
%% 2-D phase portrait
nt=length(thim);
clear thetadot
for k=1:tn
    for i=1:n-2
        for j=1:nt-1
            thetadot(j,i,k)=(theta(j+1,i,k)-theta(j,i,k))/(thim(j+1,1,k)-thim(j,1,k));
        end
    end
    
    
    figure(100)
    plot(thim(2:nt,1,k),thetadot(:,:,k))
    title(['thetadot']);
    figure(101)
    hold all
    quiver(theta(2:nt,1,k),theta(2:nt,2,k),thetadot(:,1,k),thetadot(:,2,k))
    title(['theta - thetadot direction field']);
    hold off
end

%% dp is the phase difference with ref
% d=2;
% pn=10;
dp=zeros(pn,d,tn);

figure
hold all
title(['Phase lag vs time of incoming ligt puls']);
for i=1:tn
    for j=1:d
        for k=1:pn
            dp(k,j,i)=(tim3(k,j,i)-tim3r(k,j))/tau0;
        end
        
        plot(tim3(:,j,i),dp(:,j,i))
    end
    
end
hold off
%% plot dp as tl varies: for which peak, which drop?

wp=6;%which peak?
wd=2;%which drop?


dp1=zeros(1,tn);
dp1(1,:)=dp(wp,wd,:);

% tp(1,:)=tim3(wp,1,:);%optional x-axis

figure
plot(tl,dp1)
% plot(tp,dp1)
title(['Drop ',num2str(wd),', Peak ',num2str(wp)]);



%% %%%%%%%%%%%%%%%%%%%%%%%%% old script below %%%%%%%%%%%%%%%%%%%
% i=2;
% figure
% plot(T,Y(:,5+i*9),T,Y(:,8+i*9))


%% plot [Z]

figure
hold all %hold on won't change color in loop:(
for i=1:mn
    plot(T,Y(:,3+(i)*9));
end;
hold off

%% calculate period "tau"
taumx=zeros(pn-1,mn,tn);
for k=1:tn
    for i=1:mn
        for j=1:pn-1
            taumx(j,i,k)=tim3(j+1,i,k)-tim3(j,i,k);
        end
    end
end

time3(1:pn-1,:,:)=tim3(1:pn-1,:,:);%time for plot
for k=1:tn
    figure(k)
    plot(time3(:,:,k),taumx(:,:,k));
    title(['Period vs. time, tao=',num2str(taop(k)),'s.']);
end

%% Calculate phase shift "phi"
%it's better to define the size of a matrix first
simdphi=zeros(pn-1,mn-1,tn);
simnphi=zeros(pn-1,mn-1,tn);
for k=1:tn
    for i=1:mn-1
        for j=1:pn-1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 if i<mn %locs is not the right time, tim is.
            %                      %%%%%%%%%mod 0.5
            %                      %                         simphimx(j,i)=min(abs(tim3(j,i+1)-tim3(j,i))/taumx(j,i), abs(1-abs(tim3(j,i+1)-tim3(j,i))/taumx(j,i)));
            %                      %                    else simphimx(j,i)=min(abs(tim3(j,i)-tim3(j,1))/taumx(j,i),abs(1-abs(tim3(j,i)-tim3(j,1))/taumx(j,i)));
            %                      %%%%%%%%%mod 1
            %                     simphimx(j,i)=min(abs(tim3(j,i+1)-tim3(j,i))/taumx(j,i), abs(2-abs(tim3(j,i+1)-tim3(j,i))/taumx(j,i)));
            %                 else simphimx(j,i)=min(abs(tim3(j,i)-tim3(j,1))/taumx(j,i),abs(2-abs(tim3(j,i)-tim3(j,1))/taumx(j,i)));
            %                 end;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%drop 1 as ref
            %         simdphi(j,i)=min(abs(tim3(j,i)-tim3(j,1))/taumx(j,i),abs(2-abs(tim3(j,i)-tim3(j,1))/taumx(j,i)));
            %%%%%%%%%to avoid ambiguity when rebuilding stplot from dphi
            %%%%%%%%%%%LOCAL phase difference between neighboring drops
            nphi=(tim3(j,i+1,k)-tim3(j,i,k))/taumx(j,i,k);
            while nphi<0 || nphi>1
                if nphi<0
                    nphi=nphi+1;
                else
                    nphi=nphi-1;
                end
            end
            simnphi(j,i,k)=nphi;
            %%%%%%%%%GLOBAL phase difference between drops to ref drop1
            dphi=round(tim3(j,i+1,k)-tim3(j,1,k))/taumx(j,i+1,k);
            while dphi<0 || dphi>1
                if dphi<0
                    dphi=dphi+1;
                else
                    dphi=dphi-1;
                end
            end
            simdphi(j,i,k)=dphi;
        end
        
        %         figure(tn+k)
        %         plot(time3(:,2:mn,k),simnphi(:,:,k))
        %         title(['Pair phase shift vs. time, tao=',num2str(taop(k)),'s.']);
        figure(2*tn+k)
        plot(time3(:,2:mn,k),simdphi(:,:,k))
        title(['Ref phase shift vs. time, tao=',num2str(taop(k)),'s.']);
    end
end


%% a "space-time" plot - shifted concentrations
[rowY,colY]=size(Y);
Z=zeros(rowY,mn);
figure
hold all
for i=1:mn
    Z(:,i)=Y(:,3+9*(i))+(i-1)*(3e-3);
end;
plot(T,Z);
title({'"Space-time" plot'});
hold off

%% generate a space-time plot (black bars)
% for i=1:mn
%     Fe(:,i)=Y(:,3+(i)*9);
% end;
% imagesc(mn,T,Fe) % quick and dirty
tin3 = round(tim3);%round time to integer

pnb=1;%peak number beginning, don't start from 4int
pns=pn-1;%peak number spand. don't use up all peaks
lt=5;%line thickness(half),too thick could be trouble--no larger than first peak
sec1=min([tin3(pnb,:)])-lt;
sec2=max([tin3(pnb+pns,:)])+lt;
wd=50;%width of the drop

clear stp
stp=ones(sec2-sec1+1,mn*wd);
for i=1:mn
    for j=pnb:pnb+pns
        tp=tin3(j,i);%time on the peak
        stp(tp-sec1+1-lt:tp-sec1+1+lt,(i-1)*wd+1:i*wd)=0;%more robust timewise
        if rem(j,4)==0;%&&rem(i+1,2)==0; %here comes the marker trick-for same # of oscillation
            stp(tp-sec1+1-0.5*wd:tp-sec1+1+0.5*wd,(i-1)*wd+0.5*wd-lt:(i-1)*wd+0.5*wd+lt)=0;
        end;
    end;
end;

img=mat2gray(stp,[0 1]);
imwrite(img,'a200b50ma30mM.png','png')%it seems png is tougher to zoom and shrink
%% phase shift as a function of space
figure
plot(simnphi(pn-1,:))
%% make an animation
% pn=230;
% n=27;
for j=1:pn-1
    figure(1); clf; %Clear current figure window
    set(gcf,'DoubleBuffer','on'); % somehow this works to avoid annoying flashes when animating graphs
    axis([0 n 0 1]);
    
    title(['# of oscillations: ',num2str(j)],'Color','b')
    hold on
    plot(simnphi(j,:),'Marker','o','LineWidth',1,'Color','b');
    %     plot(simnphia100b0(j,:),'Marker','o','LineWidth',1,'Color','b');
    %     plot(simnphia100b100(j,:),'Marker','s','LineWidth',1,'Color','g');
    %     plot(simnphia100b500(j,:),'Marker','x','LineWidth',1,'Color','r');
    %     plot(simphimx2000(j,:),'Marker','+','LineWidth',1,'Color','black');
    M(j) = getframe;
end;
hold off
movie(M)
movie2avi(M, 'a100b0-500ma640.avi', 'compression', 'None');

%% phase shift as a function of space (odd even separate)
for i=1:19
    %drop1 as ref, plot start from drop2
    edphi(1,i)=simnphi(pn-1,2*i-1);%drop2,4,...
    odphi(1,i)=simnphi(pn-1,2*i);%drop3,5,...
end

npair=1:19;

figure
hold all
plot(2*npair-1,edphi)
plot(2*npair,odphi)
plot(simnphi(pn-1,:))
hold off

%% Reaction or diffusion? du/dt
h = 0.16; %[H+](Mole)
A = 0.3; %[BrO3-]
m = 0.4; %[MA]
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

lw = 0.02; %length water drop(cm)~200um
lo = 0.005; %length oil drop(cm)~50um
D = 1e-5; %diffusion coefficient cm^2/s, I assume it's the same for both Br2 and radical
pB = 2.5;% partition coefficient for [Br2], Co*/Cw* at equilibrium,
pR = 1;% partition coefficient for [BrO2*], NOT sure about the value:(
sl = lw + lo;% sum of water drop length lw and oil drop length lo
cdw = 2*D/(lw*sl);%effective reaction rate due to diffusion, water to oil
cdo = 2*D/(lo*sl);%effective reaction rate due to diffusion, oil to water
bc=0.05;

% tl=zeros(1,n);%time of light ptb
%%bdy drop
% tl(1)=1e6;
% tl(n)=1e6;
% %%light perturbation
%
% tl(2)=110;
% tl(4)=110;

% lp=1e-4;%light on %1e-4 for high ma, 10e-4 for low ma
% % lp=0;%light off
% for i=1:n
%     for j=1:51507
%         ki(j,i)=lp*(1-sign(T(j,1)-tl(i)))*0.5;
%     end;
% end;

i=1;
% dcdt(:,1) = -kred*Y(:,6+i*9).*Y(:,7+i*9)+k9*Y(:,3+i*9)*m+k10*Y(:,3+i*9)*m-ki(j,i+1).*Y(:,7+i*9)*b/(bc+b);
dudt(:,1)= k5*Y(:,2+i*9).*Y(:,4+i*9)*h-k6*Y(:,5+i*9)-k7*Y(:,5+i*9)*m + cdw*(Y(:,8+(i-1)*9)+Y(:,8+i*9)-2*pB*Y(:,5+i*9));
dudtR(:,1)= k5*Y(:,2+i*9).*Y(:,4+i*9)*h-k6*Y(:,5+i*9)-k7*Y(:,5+i*9)*m;
dudtD(:,1)= cdw*(Y(:,8+(i-1)*9)+Y(:,8+i*9)-2*pB*Y(:,5+i*9));

% td=size(T);
% for j=1:td(1)
%     dudt(j,1)=k5*Y(j,2+i*9)*Y(j,4+i*9)*h-k6*Y(j,5+i*9)-k7*Y(j,5+i*9)*m + cdw*(Y(j,8+(i-1)*9)+Y(j,8+i*9)-2*pB*Y(j,5+i*9));
%     dudtR(j,1)=k5*Y(j,2+i*9)*Y(j,4+i*9)*h-k6*Y(j,5+i*9)-k7*Y(j,5+i*9)*m;
% end
figure
plot(T,dudt,T,dudtR,T,dudtD)
% figure
% plot(T,dcdt)
%% [Br2] in BZ and oil
i=2;
figure
plot(T,Y(:,5+i*9),T,Y(:,8+i*9))