% start simulation near steady state (ss)
% and compare with LSA result at early time and small ptb
close all
clear




rng('shuffle') % IMPORTANT:seeds the random number generator based on the current time
global k1 k2 k3 k4 k7 k9 k10 c0 cmin

abstol=1e-13;%1e-9;%


h = 0.16; %[H+](Mole)
A = 0.3; %[BrO3-]
n = 40;

ncase=6;
dnm=2;%drop displayed
eps=1e-3;%1e-1;%
cutoff=1e-2;

if ncase==1
    qmax=0;
    m=0.001;%case a (lsa)
% m=0.065;%mu=29;
    %case a
%     eigvec=real([0.121919885509685;-0.591222614746859;0.693466853992190;0.393312959033092;]);
        eigvec=real([-0.0836596049016332;-0.791241636857404;0.295129188761342;0.528995750977256;]);%oxidized--weird
    %     warning: Failure at t=4.058768e-003.  Unable to meet
    % integcutoffn tolerances without reducing the step
    % size below the smallest value allowed
    % (1.387779e-017) at time t.
%     eigvec=real([0.114838426145575;-0.432559795743582;0.884164318913889;0.134006030260644;]);%m=0.065;mu=29;
    
elseif ncase==2
    qmax=0;
    m=0.9;%case b
    %case b
    eigvec=real([-0.0477548056565754 + 0.0266985379682336i;0.252186709169694 - 0.235490733687297i;-0.936961558542126 + -0.00000000000000i;-0.00680429711075779 + 0.00306333605369280i;]);
%     eigvec1=[-0.0477548056565754 + 0.0266985379682336i;0.252186709169694 - 0.235490733687297i;-0.936961558542126 + -0.00000000000000i;-0.00680429711075779 + 0.00306333605369280i;];
%     eigvec2=[-0.0477548056565754 - 0.0266985379682336i;0.252186709169694 + 0.235490733687297i;-0.936961558542126 + -0.00000000000000i;-0.00680429711075779 - 0.00306333605369280i;];
elseif ncase==3
    qmax=1;
    m=0.2;%case c, long term pi-s state
%         m=0.05;
%         m=0.4;
%         m=0.4;%change mu to 0.15
    %case c
    eigvec=real([0.0548422288208566;-0.501496360621181;0.863194872476654;0.0197064038242457;]);
%         eigvec=real([-0.124968894223853;0.574534952224631;-0.805948983933662;-0.0688374857545330;]);
%         eigvec=real([-0.0501314070911387;0.410784174350276;-0.910283271021722;-0.0112947166595114;]);
%         eigvec=real([-0.0412409776046593;0.291779462324381;-0.955527896609906;-0.0114177901608172;]);%for mu=0.15
elseif ncase==4
%         qmax=0.4;
    qmax=0.35;
    m=0.005;
%         eigvec=real([-0.125412200885052;0.625613184483952;-0.746575275808960;-0.188428450133710;]);%unphysical
    eigvec=real([-0.130210016876501;0.613487668020959;-0.752050592888004;-0.202726758046460;]);%ox
elseif ncase==5
    qmax=0.5;
    m=0.75;%case e
    %case e qmax=0.5
    eigvec=real([-0.0439871271736746 + 0.0199940082773635i;0.266570466415269 - 0.180671997654321i;-0.945467954730988 + -0.00000000000000i;-0.00684706191899805 + 0.00257920426361028i;]);
    %400 drops LSA (very close result to 40 drops):
    %     eigvec=real([-0.0439860093422857 + 0.0199441899496342i;0.266897083484680 - 0.180322923180042i;-0.945443643260192 + -0.00000000000000i;-0.00683840475562073 + 0.00256984666117376i;]);
elseif ncase==6
    qmax=1;
    m=0.6;%case f
    %case f
    eigvec=real([-0.0407684391227943 + 0.00883251395375749i;0.283091337455353 - 0.0823146893528653i;-0.954616508929856 + -0.00000000000000i;-0.00700998242508097 + 0.00129948766949407i;]);
    
end

k1 = 2e6*h; %(1/Molsec)=(1/Ms)
k2 = 2*A*h.^2; %(1/s)
k3 = 3000; %(1/Ms)
k4 = 42*A*h; %(1/s)
k7 = 29*m; %(1/s)

%     k9p=0.07;%0.12;% for m>+0.1;0.07 for m<0.1;% (1/s),k9p=0~0.2

if m>0.1
    k9p=0.12;
else
    k9p=0.07;
end

k9 = k9p*m;
k10 = 0.05*m; %(1/s)
kr = 2e8; %(1/Ms)
kred = 5e6; %(1/Ms)

c0=4.2e-3;%3e-3;%
cmin=sqrt(2*kr*(k9+k10)*c0/kred^2);


%  get the symbolic Jacobian
syms y1 y2 y3 y4

f1 =[-k1*y1*y2+k2*y2-2*k3*y1^2+k4*y1*(c0-y3)/(c0-y3+cmin);
    -3*k1*y1*y2-2*k2*y2-k3*y1^2+k7*y4+k9*y3;
    2*k4*y1*(c0-y3)/(c0-y3+cmin)-k9*y3-k10*y3;
    2*k1*y1*y2+k2*y2+k3*y1^2-k7*y4];
v1=[y1, y2, y3, y4];
J=jacobian(f1,v1);

% solve for stead state X*, Y*, Z*, U*
v0=1e-4*[1 1 1 1];%trial solution
% v0=[5.92*1e-7 1.44e-6 0.00419924 0.0000195482];%case a ma=0.001;
% v0=[6.16443e-7 6.59949e-6 0.00400912 0.0000186632];%case d ma=0.005;
options=optimset('TolFun',abstol);% <or= AbsTol for ODE solver
[v,fval] = fsolve(@vefunc,v0,options);
format long %more accurate
disp('Steady State: X*, Y*, Z*, U*=')
disp(v)
%eigenvalue at steady state
y1=v(1); y2=v(2); y3=v(3); y4=v(4);%@steady state
J1=subs(J);%replace all symbols with numbers
eig(J1);%negative value means unstable

%% IC=SS+pt
runtime=3000;%sec
timestep=1e-1;%1e-3;% %sec
% eps=1e-3;%1e-5;%

dc=zeros(4,1);

dc(4)=0.6;%29;%0.15;%
dc(1)=dc(4)/50;
%                 bz=round(1e4*sqrt((1/dc(4))*1e-5*2.5/2));
%                 oil=bz;



yss=[y1;y2;y3;y4];
nyss=norm(yss);

% norm(eigenvector)=1
ic=zeros(1,4*n);
for i=1:n
    ic(1+4*(i-1):4+4*(i-1))=yss+eps*nyss*cos(i*qmax*pi)*eigvec;
end

%random IC SS20, double checked with phase plot
% cvy=20/100;
% icy=zeros(4,n);
% for i=1:4
%     icy(i,:)=abs(normrnd(yss(i),yss(i)*cvy,1,n));
% end
% 
% ic=zeros(1,4*n);
% for i=1:n
%     ic(1+4*(i-1):4+4*(i-1))=icy(1:4,i);%1);%
% end

% dnm=4;

clear options
options = odeset('RelTol',abstol*10,'AbsTol',abstol,'OutputFcn', @odephas3, 'OutputSel',[1 2 3]+(dnm-1)*4);

clear T Y

[T,Y] = ode15s(@(t,y) fourvar(t,y,n,dc,h,A,m),[0:timestep:runtime],ic,options);
nT=length(T);

%% realm of LSA validity
% look at one drop on the ring
theta=zeros(5,nT);
Yss=transpose(yss);

for i=1:nT
    for dn=1:5
        theta(dn,i)=acos(dot(Y(i,1+(dn-1)*4:4+(dn-1)*4)-Yss,eps*nyss*eigvec*cos(dn*qmax*pi))/norm(Y(i,1+(dn-1)*4:4+(dn-1)*4)-Yss)/norm(eps*nyss*eigvec*cos(dn*qmax*pi)));
    end
end

theta=theta*180/pi;%radian to degree
figure
subplot(2,1,1)
plot(T,theta)
ylabel('\theta')
subplot(2,1,2)
plot(T,Y(:,1+(dnm-1)*4:4+(dnm-1)*4))
ylabel(['conc. in drop ', num2str(dnm)])
xlabel('time (s)')

%% when it hits limit cycle? 1
% Y(t)-Y(t+tau) approach zero
clear pks locs
[pks,locs]=findpeaks(Y(:,7),'minpeakheight',1e-3,'minpeakdistance',100);
tau=locs(end)-locs(end-1);%in number of time steps, not in seconds
Ydn=Y(:,1+(dnm-1)*4:4+(dnm-1)*4);%for dnm drop
mnlc=zeros(1,nT);
% cutoff=5e-3;
for i=1:nT-tau
    mnlc(i)=norm(Ydn(i,:)-Ydn(i+tau,:));
end
figure
plot(T,mnlc)
title('difference with one lc period later')
ylabel('|Y(t)-Y(t+tau)|')
xlabel('time(s)')
% threshold: th= mnlc(t)/mnlc(t=0)=1%
t90=find(theta(dnm,:)>90,1);
t901=t90*timestep% # of time step
th=find(mnlc<mnlc(1)*cutoff,1);
th1=th*timestep
lv1=t90/th

%% when it hits limit cycle? 2
% Y(t)-Ylc approach zero
clear pks locs
[pks,locs]=findpeaks(Y(:,7),'minpeakheight',1e-3,'minpeakdistance',100);
tau=locs(end)-locs(end-1);%in number of time steps, not in seconds
Ylc=Y(locs(end-1):locs(end),1+(dnm-1)*4:4+(dnm-1)*4);
Ydn=Y(:,1+(dnm-1)*4:4+(dnm-1)*4);%for dnm drop
mnlc=zeros(1,nT);%minimum norm from limit cycle
nlc=zeros(1,tau+1);% norm from limit cycle

for j=1:nT
    for i=1:tau+1
        nlc(i)=norm(Ydn(j,:)-Ylc(i,:));
    end
    mnlc(j)=min(nlc);
end

figure
plot(T,mnlc)
title('minimum distance to limit cycle')
xlabel('time(s)')
ylabel('minimum norm to limit cycle')
% threshold: th= mnlc(t)/mnlc(t=0)=1%
t90=find(theta(dnm,:)>90,1);% # of time step
t902=t90*timestep
th=find(mnlc<mnlc(1)*cutoff,1);
th2=th*timestep
lv2=t90/th

lveps=norm(Ydn(t90,:)-Yss)/norm(Yss)%the largest perturbation within linear validity

%% FFT
Yf=zeros(n,nT,4);


for j=1:nT
    for i=1:n
        for k=1:4
            Yf(i,j,k)=Y(j,k+(i-1)*4);
        end
    end
end


ewqt=zeros(nT,4);
% eiqx=zeros(1,n);
% for i=1:n
%     eiqx(i)=cos(i*qmax*pi);
% end
delct=zeros(nT,4);
delc0=zeros(1,4);
for j=1:nT
    for k=1:4
       
        if qmax==0
            delct(j,k)=Yf(1,j,k)-yss(k);
            delc0(k)=Yf(1,1,k)-yss(k);
        elseif qmax==1
            delct(j,k)=Yf(2,j,k)-yss(k);
            delc0(k)=Yf(2,1,k)-yss(k);
        elseif isinteger(2/qmax)==1
            delct(j,k)=Yf(2/qmax,j,k)-yss(k);%so that cos(i*qm*pi)=1, but it's picky about qmax
            delc0(k)=Yf(2/qmax,1,k)-yss(k);
        else
            delct(j,k)=Yf(n,j,k)-yss(k);%so that cos(i*qm*pi)=1
            delc0(k)=Yf(n,1,k)-yss(k);
        end
        ewqt(j,k)=delct(j,k)/delc0(k);
    end
end
%%%%%%%%%%for qmax%%%%%%%%%%%
wqmt=zeros(nT,4);
for k=1:4
    wqmt(:,k)=log(ewqt(:,k));
end

wqm=zeros(nT-1,4);
for k=1:4
    for j=1:nT-1
        wqm(j,k)=(wqmt(j+1,k)-wqmt(j,k))/timestep;%slope,d/dt
    end
end
%%%%%%%%%%for all q%%%%%%%%%%%
% wqt=log(ewqt);
% wq=zeros(n,nT-1,4);
% for k=1:4
%     for i=1:n
%         for j=1:nT-1
%             wq(i,j,k)=(wqt(i,j+1,k)-wqt(i,j,k))/timestep;%slope,d/dt
%         end
%     end
% end

%% for qmax
figure
plot(T,wqmt)
xlabel('t(s)')
ylabel('w_qt')
title(['MA=',num2str(m),'M, \mu_U=',num2str(dc(4))])

figure
plot(T(1:end-1),wqm)
xlabel('t(s)')
ylabel('w_q')
title(['MA=',num2str(m),'M, \mu_U=',num2str(dc(4))])


%% for all q
% figure(1)
% qmax=1;
% wqk=zeros(n,nT-1);
% for k=1:4
%     wqk=wq(:,:,k);
%     hold all
%     plot(T(1:end-1),wqk(qmax*n/2+1,:))
%     xlabel('t(s)')
%     ylabel('w_q')
%     title(['MA=',num2str(m),'M, \mu_U=',num2str(dc(4))])
% end
% hold off
%
% figure(2)
% wqtk=zeros(n,nT);
% for k=1:4
%     wqtk=wqt(:,:,k);
%     hold all
%     plot(T(1:end),wqtk(qmax*n/2+1,:))
%     xlabel('t(s)')
%     ylabel('w_qt')
%     title(['MA=',num2str(m),'M, \mu_U=',num2str(dc(4))])
% end
% hold off
%
% figure(3)
% q=(0:20)*2/n;
% for k=1:4
%     subplot(2,2,k)
%     plot(q,wq(1:21,1:(runtime/10/timestep):end,k))
%     xlabel('q')
%     ylabel('w_q')
%     title(['MA=',num2str(m),'M, \mu_U=',num2str(dc(4)),', C_',num2str(k)])
% end
%% get simulated eigenvector at t
%case a
% sig=real([0.708877268611527;]);
%case b
% sig=real([0.0476489869309855 - 0.112176244123156i;]);
%case c
% sig=real([0.219260965511917;]);
% eigv=zeros(4,1);
% for k=1:4
%     eigv(k)=(Yf(1,end,k)-yss(k))/eps/nyss;
%     %     eigv(k)=(Yf(1,end,k)-yss(k))/eps/yss(k);
% end
% disp('relative difference from LSA eigenvector:')
% disp((eigv-eigvec)./eigvec)
%% Long term profile
% clear options
% abstol=1e-9;
% options = odeset('RelTol',abstol*10,'AbsTol',abstol,'OutputFcn', @odephas3, 'OutputSel',[1 2 3]);
% runtime1=500;
% timestep1=1e-1;
% % clear ic
% % ic=1e-4*rand(1,4*n);
% clear T1 Y1
% [T1,Y1] = ode15s(@(t,y) fourvar(t,y,n,dc,h,A,m),[0:timestep1:runtime1],ic,options);
% % odephas3 ?Three-dimensional phase plane plotting
% % [1 2 3] means y(1) y(2) y(3)
% Y2=zeros(length(T1),n);
% for i=1:n
%     Y2(:,i)=Y1(:,3+(i-1)*4)+c0*(i-1);
% end
% figure(4)
% plot(T1,Y2)
% xlabel('t(s)')
%% plot concentcutoffn - SS
% Yss=zeros(nT,4);
% 
% figure(5)
% for k=1:4
%     subplot(2,2,k)
%     plot(T,delct(:,k));
%     title(['\deltaC_',num2str(k)]);
%     xlabel('time(s)')
% end

%% 2D limit cycle
% close all
% 
% x=Y1(:,1);
% y=Y1(:,2);
% z=Y1(:,3);
% % u=Y1(:,4);
% 
% figure(1)
% plot3(x,y,z)
% xlabel('x')
% ylabel('y')
% zlabel('z')

% figure(1)
% plot(x,y)
% xlabel('x')
% ylabel('y')
% figure(2)
% plot(x,z)
% xlabel('x')
% ylabel('z')
% figure(3)
% plot(x,u)
% xlabel('x')
% ylabel('u')
% figure(4)
% plot(y,z)
% xlabel('y')
% ylabel('z')
% figure(5)
% plot(y,u)
% xlabel('y')
% ylabel('u')
% figure(6)
% plot(z,u)
% xlabel('z')
% ylabel('u')

%% principal component analysis
% close all
% [coeff,score,latent]=princomp(Y1(:,1:4));

%The columns of coeff are in order of decreasing component variance.
%[COEFF,SCORE] = princomp(X) returns SCORE, the principal component scores; that is, the representation of X in the principal component space. Rows of SCORE correspond to observations, columns to components.
% cumsum(latent)./sum(latent)
%This show that two components account for 99.9% of the variance:

% x=score(:,1);
% y=score(:,2);
% z=score(:,3);
% u=score(:,4);

% 
% figure(2)
% plot3(x,y,z)
% xlabel('1st component')
% ylabel('2nd component')
% zlabel('3rd component')

% figure(2)
% plot(x,z)
% xlabel('x')
% ylabel('z')
% figure(3)
% plot(x,u)
% xlabel('x')
% ylabel('u')
% figure(4)
% plot(y,z)
% xlabel('y')
% ylabel('z')
% figure(5)
% plot(y,u)
% xlabel('y')
% ylabel('u')
% figure(6)
% plot(z,u)
% xlabel('z')
% ylabel('u')
%% double check
% sigt= log(abs(fft(Y2(end,:)-yss(2)))./abs(fft(Y2(1,:)-yss(2))));
% sigt(1)-wqt(1,end,2)%should give zero

%% load comsol
% comsol=load('20drops_C3_pX001.txt');
% dn=20;
% tn=length(comsol)/dn;
% c0=4.2e-3;
% comT=0:tn-1;
% for i=1:dn
%     comC3(1:tn,i)=comsol(1+(i-1)*tn:tn+(i-1)*tn,2)+c0*(i-1);
% end
% figure
% plot(comT,comC3);