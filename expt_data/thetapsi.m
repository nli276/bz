% partyly taken from phasematrix.m for calculation of theta and psi
%% input
d=3; % # of droplets

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
            psimx(j,i)=abs(theta(j,i+1)-theta(j,i));
            
        else psimx(j,i)=abs(theta(j,1)-theta(j,i));
        end;
    end;
end;

psimx=psimx/(2*pi);