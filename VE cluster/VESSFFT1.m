% start simulation near steady state (ss)
% FFT of Turing pattern as a function of time
% close all
% clear

function [qmx,lmx,wmx,qmxn,lmxn,wmxn]=VESSFFT1(myseed,mytaskid)
tic;
rng(myseed)
% rng('shuffle') % IMPORTANT:seeds the random number generator based on the current time
global k1 k2 k3 k4 k7 k9 k10 c0 cmin
nsz=30;
nma=46;

qmx=zeros(nsz,nma);
lmx=zeros(nsz,nma);%wavelength
wmx=zeros(nsz,nma);
qmxn=zeros(nsz,nma);%n means nonlinear
lmxn=zeros(nsz,nma);
wmxn=zeros(nsz,nma);
h = 0.16; %[H+](Mole)
A = 0.3; %[BrO3-]
for malp=1:nma
%     if malp<10
%         m=0.001*malp;
%     elseif malp<20
%         m=0.01+0.01*(malp-10);
%     else
%         m=0.2+0.1*(malp-20);
%     end
    m=0.001+0.02*(malp-1);
    
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
    v0=[1e-4 1e-4 1e-4 1e-4];%trial solution
    options=optimset('TolFun',1e-9);% <or= AbsTol for ODE solver
    [v,fval] = fsolve(@vefunc,v0,options);
    format long %more accurate
    %     disp('Steady State: X*, Y*, Z*, U*=')
    %     disp(v)
    %eigenvalue at steady state
    y1=v(1); y2=v(2); y3=v(3); y4=v(4);%@steady state
    J1=subs(J);%replace all symbols with numbers
    eig(J1);%negative value means unstable
    
    %% SS as IC
    runtime=6000;%10;%
    timestep=1;%1e-3;%
    n=40;%100;%# of droplets on a ring, n>=2
    %     nc=1;%10;%# of cycles
    
    for szlp=1:nsz
        %%%%%%%%%%%%%%%geometry%%%%%%%%%%%%
        %         bz = 0.001*szlp; %length of bz (water) drop(cm), 0.02 is 200um
        %         oil = bz;%0.015; %length of oil drop(cm)
        
        
        %%%%%%%%%%%%diffusion%%%%%%%%%%%%%%%
        %         par=[0.05 0 0 2.5];%[0 0 0 0];%[Px Py Pz Pu]
        %zero for ref partition coefficient for the 4 species, i.e. drops are not coupled for ref
                
%         cd=zeros(4,1);
        dc=zeros(4,1);
        %         D = 1e-5; %diffusion coefficient cm^2/s, I assume it's the same for all species
        %         for i=1:4
        %             dc(i)=D*par(i)/(bz*(bz+oil));
        %         end
%      
%         cd(4)=0.02+szlp-1;
%         cd(1)=cd(4)*50;
        %we don't want dc(2) and dc(3) to be Inf
%         dc(4)=1/(0.02+szlp-1);
        dc(4)=0.02*szlp;
        dc(1)=dc(4)/50;
%         bz=round(1e4*sqrt((0.02+szlp-1)*1e-5*2.5/2));
%         oil=bz;

        % check if ss is really ss - checked
        % ic=zeros(1,4*n);
        % for i=1:n
        %     ic(1+(i-1)*4)=v(1);
        %     ic(2+(i-1)*4)=v(2);
        %     ic(3+(i-1)*4)=v(3);
        %     ic(4+(i-1)*4)=v(4);
        % end
        
        
        
        % SS + 20% fluctuation as IC
        cvy=20/100; %The coefficient of variation (CV) equals the standard deviation divided by the mean (expressed as a percent).
        yss=[y1;y2;y3;y4];
        icy=zeros(4,n);
        for i=1:4
            icy(i,:)=abs(normrnd(yss(i),yss(i)*cvy,1,n));
        end
        
        ic=zeros(1,4*n);
        for i=1:n
            ic(1+4*(i-1):4+4*(i-1))=icy(1:4,i);%1);%
        end
        
        clear options
        options = odeset('RelTol',1e-8,'AbsTol',1e-9);
        
        clear T Y
        
        [T,Y] = ode15s(@(t,y) fourvar(t,y,n,dc,h,A,m),[0:timestep:runtime],ic,options);
        
        %generate space time plot png image
        s=3;
%         Fe=zeros(runtime+1,n*40);%30 pixel for drop, 10 pixel for gap (black)
%         for i=1:n
%             for j=1:runtime+1
%                 Fe(j,(i-1)*40+1:(i-1)*40+30)=Y(1+(j-1)/timestep,s+(i-1)*4);
%             end
%         end
        
%         img=mat2gray(transpose(Fe),[0 c0]);
%         filename=['a', num2str(bz), 'b',num2str(oil), 'umMA',num2str(m*1e3),'mM',num2str(n),'drops','.png'];
        %             imwrite(img,filename);
        %             figure(1)
        %             imshow(filename)
        %
        %% FFT
        Yf=zeros(n,length(T));%for FFT
        %                         Yh=zeros(n,length(T));%for histrogram
        
        %binarize Z(n)
        
        %                         for j=1:length(T)
        %                             for i=1:n
        %                                 if Y(j,s+(i-1)*4)>=0.99*c0
        %                                     Yh(i,j)=1;
        %                                 end
        %                             end
        %                         end
        
        %fft for Z(n)
        for j=1:length(T)
            for i=1:n
                if Y(j,s+(i-1)*4)>0.999*c0
                    Yf(i,j)=c0;
                else
                    Yf(i,j)=Y(j,s+(i-1)*4);%transpose for FFT(t)
                end
            end
        end
        Yf=Yf/c0;%normalize
  
        Yln=Yf(:,1:10);%linear
        Ynl=Yf(:,end-5000:end);%nonlinear
        
        
        %
        ft2d=fft2(Yln);
        ft2d(1,1)=0;
        qs=1;%1 or 2
        ws=1;%1 or 2
        absft2d=abs(ft2d(qs:n/2+1,:));%(1,1) is at maximum always, but this maximun isn't always meaningful
                
        [C1,I1]=max(absft2d);
        [C2,I2]=max(C1);
        qmx(szlp,malp)=(I1(I2)-2+qs)/(n/2);
        lmx(szlp,malp)=n/(I1(I2)-2+qs);%wavelength in drop number
        wmx(szlp,malp)=I2-2+ws;
        
        ft2dn=fft2(Ynl);
        ft2dn(1,1)=0;
        qs=1;%1 or 2
        ws=1;%1 or 2
        absft2dn=abs(ft2dn(qs:n/2+1,ws:30));%(1,1) is at maximum always, but this maximun isn't always meaningful
                
        [C3,I3]=max(absft2dn);
        [C4,I4]=max(C3);
        qmxn(szlp,malp)=(I3(I4)-2+qs)/(n/2);
        lmxn(szlp,malp)=n/(I3(I4)-2+qs);%wavelength in drop number
        wmxn(szlp,malp)=I4-2+ws;
        
        %         h1=figure(1);
        %         subplot(2,1,1)
        %         imagesc(img)
        %         xlabel('t')
        %         ylabel('d')
        %         title(['a=', num2str(bz*1e4), ',b=',num2str(oil*1e4), 'um, MA=',num2str(m*1e3),'mM,',num2str(n),'drops']);
        %
        %         subplot(2,1,2)
        %         imagesc(absft2d)
        %         xlabel('w')
        %         ylabel('q')
        %         title(['q=', num2str((I1(I2)-2+qs)/(n/2)),', l=',num2str(n/(I1(I2)-2+qs)),', w=',num2str(I2-2+ws)]);
        %         saveas(h1,filename,'png')
        %       pause
        %       end
        
    end
end
matname=['qlwmx_',num2str(mytaskid),'.mat'];
save(matname,'qmx','lmx','wmx','qmxn','lmxn','wmxn')
% qmx
% lmx
% wmx
toc;
