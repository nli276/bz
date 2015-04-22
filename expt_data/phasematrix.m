%Written by Ning Li on Sep 29, 2009 in order to
%calculate phase (theta) for each droplet and do correlation stuff.


%% input
d=3; % # of droplets
%f=8231; % # of frames


%% cut data to same length
for i=1:d
    eval(['n(i)=length(locsy' num2str(i) ');']);
    nmin=min(n);
end;

%% put locsy# variable into a matrix
for i=1:d
    eval(['locsmx(1:nmin,i)=locsy' num2str(i) '(1:nmin);']);
end;

%% the frame # of the first end peak of all droplets, 
%i.e. cut the space-time data to same (minimal) length.
minplocs=min(locsmx(nmin,:));
%display (minplocs)

%% calculate theta, each column represent a droplet
for i=1:d
    j=1;
    for k=1:nmin-1
        
        while j<minplocs
            
            if locsmx(k,i)<=j-1 && j-1<locsmx(k+1,i)
                theta(j,i)=2*pi*(j-1-locsmx(k,i))/(locsmx(k+1,i)-locsmx(k,i));
            elseif j-1<locsmx(k,i);
                theta(j,i)=0;
            else break
            end;
            
            j=j+1;
            
        end;
        
    end;
    
end;

%% calculate theta, each column represent a droplet
for i=1:d
    j=1;
    for k=1:nmin-1
        
        while j<minplocs
            
            if locsmx(k,i)<=j-1 && j-1<locsmx(k+1,i)
                theta(j,i)=2*pi*(j-1-locsmx(k,i))/(locsmx(k+1,i)-locsmx(k,i))+2*pi*(k-1);
            elseif j-1<locsmx(k,i);
                theta(j,i)=0;
            else break
            end;
            
            j=j+1;
            
        end;
        
    end;
    
end;
%% Calculate phase shift psi using theta (peak locs -> phi)
for i=1:d
    for j=1:minplocs-1
        if i<d
            psimx(j,i)=min(abs(theta(j,i+1)-theta(j,i)),abs(2*pi-abs(theta(j,i+1)-theta(j,i))));
            
        else psimx(j,i)=min(abs(theta(j,1)-theta(j,i)),abs(2*pi-abs(theta(j,1)-theta(j,i))));
        end;
    end;
end;

psimx=psimx/(2*pi);

%% Calculate phase shift psi using theta (peak locs -> phi)
for i=1:d
    for j=1:minplocs-1
        if i<d
            psimx(j,i)=abs(theta(j,i+1)-theta(j,i));
            
        else psimx(j,i)=abs(theta(j,1)-theta(j,i));
        end;
    end;
end;

psimx=psimx/(2*pi);

%% calculate spatial phase correlation (spc) as a function of time t
%and correlation length l -- C(t,l).
for j=1:minplocs-1
    for l=1:d-1
        for i=1:d-l
            sindtheta(j,i,l)=sin(theta(j,i)-theta(j,i+l));%it's L not 1
            cosdtheta(j,i,l)=cos(theta(j,i)-theta(j,i+l));
        end;
        sumsin(j,l)=sum(sindtheta(j,:,l));
        sumcos(j,l)=sum(cosdtheta(j,:,l));
        
        spc(j,l)=sqrt((sumsin(j,l)/(d-l))^2+(sumcos(j,l)/(d-l))^2);
    end;
end;

%% make a cartoon for drops in phase space (running on a circle)
clear i %it's not always sqrt(-1)...
for k=1:d
    for j=1:minplocs-1
        phase(j,k)=exp(i*theta(j,k));% use theta to make a complex number so we can plot in circle.
    end;
end;

for j=1:minplocs-1
    figure(1); clf; %Clear current figure window
    set(gcf,'DoubleBuffer','on'); % somehow this works to avoid annoying flashes when animating graphs
    axis([-1.5 +1.5 -1.5 1.5]);
    
    title(['frame # : ',num2str(j)],'Color','b')
    hold on
    plot(phase(j,:),'o');
    
end;
hold off

%% make a 2-D density plot
bin=0:pi/10:2*pi;

for i=1:100 %minplocs-1
h(i,:)=hist(theta(i,:),bin);
end;

img=mat2gray(h,[0,d]);
imwrite(img,'densityplot.bmp','bmp')