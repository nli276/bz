% linear stability analysis
%% jacobian for the VE model with diffusive terms
close all
clear


rng('shuffle') % IMPORTANT:seeds the random number generator based on the current time
global k1 k2 k3 k4 k7 k9 k10 c0 cmin

abstol=1e-13;%1e-9;%



h = 0.16; %[H+](Mole)
A = 0.3; %[BrO3-]

n=40;

% m=0.001;%case a
% m=0.065;%case a with mu=29
% m=0.002;%case d
m=0.9;%case b
% m=0.2;%case c
% m=0.005;%case d
% m=0.75;%case e
% m=0.6;%case f
% m=0.05;
% m=0.4;

k1 = 2e6*h; %(1/Molsec)=(1/Ms)
k2 = 2*A*h.^2; %(1/s)
k3 = 3000; %(1/Ms)
k4 = 42*A*h; %(1/s)
k7 = 29*m; %(1/s)

% k9p=0.12;%0.07;% for m>+0.1;0.07 for m<0.1;% (1/s),k9p=0~0.2

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
%trial solution
v0=1e-4*ones(1,4);
% v0=1e-4*rand(1,4);
% v0=c0*ones(1,4);
%case a ma=1mM, mu=0.6
% v0=[5.92*1e-7 1.44e-6 0.00419924 0.0000195482];
%case d ma=2mM mu=0.6
% v0=[6.06769e-7 2.80445e-6 0.00419719 0.0000195386];
%case b ma=0.9M 
% v0=[2.81142e-7 0.00441081 0.00430157 0.0000330034];
%case c ma=0.2
% v0=[2.81138e-7 0.000968025 0.00424818 0.0000325938];
%case d ma=0.005
% v0=[6.16443e-7 6.59949e-6 0.00400912 0.0000186632];
%case e ma=0.75
% v0=[2.81141e-7 0.00366814 0.00429276 0.0000329358];
%case f ma=0.6
% v0=[2.81141e-7 0.00292785 0.00428301 0.000032861];


options=optimset('TolFun',abstol);% <or= AbsTol for ODE solver
[v,fval] = fsolve(@vefunc,v0,options);
format long %more accurate
disp('Steady State: X*, Y*, Z*, U*=')
disp(v)
%eigenvalue at steady state
y1=v(1); y2=v(2); y3=v(3); y4=v(4);%@steady state
J1=subs(J);%replace all symbols with numbers
% eig(J1)%negative value means unstable

%% Dispersion relation


dc=zeros(4,1);

dc(4)=0.6;%28.58;%
dc(1)=dc(4)/50;
bz=round(1e4*sqrt((1/dc(4))*1e-5*2.5/2));
oil=bz;

lam=zeros(4,n/2+1);
% lams=zeros(1,n/2+1);
q=zeros(1,n/2+1);
eigvect=zeros(4,n/2+1);
eigval=zeros(1,n/2+1);
for j=1:n/2+1
    syms lambda
    q(j)=(2/n)*(j-1);
    diffmx=[4*dc(1)*sin(q(j)*pi/2)^2 0 0 0;0 0 0 0;0 0 0 0;0 0 0 4*dc(4)*sin(q(j)*pi/2)^2];
    JD=J1-diffmx;
    [vjd,djd]=eigs(JD);
    clear I
    [~,I]=max(max(real(djd)));
    eigval(j)=djd(I,I);
    eigvect(:,j)=vjd(:,I);
    detJD=det(JD-lambda*eye(4,4));
    lambda=solve(detJD);
    lam(:,j)=subs(lambda);
    
%     clear J
%     [~,J]=max(real(lam(:,j)));
%     lams(1,j)=lam(J,j);
    
    clear lambda
end

lamsort=sort(lam);
lamRe=real(lamsort);
lamIm=imag(lamsort);

figure(1)
for i=1:4
subplot(2,4,i)
plot(q,lamRe(i,:))
ylabel('Re\lambda')
xlabel('q')
subplot(2,4,i+4)
plot(q,lamIm(i,:))
ylabel('Im\lambda')
xlabel('q')
end

% figure(2)
% plot(q,lams)
%the problem Im(lams) looks messy is due to two equal Re(lam) sometimes...

% lamodu=abs(lamsort);
% figure(3)
% plot(q,lamodu)
% xlabel('q')
% ylabel('|\lambda|')
