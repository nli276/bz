%**created by Ning Li on Feb 4, 2010 for calculate periods of oscillations and phase shift between droplets
t = 5; %total number of profiles (droplets)

%% cut the data to same length
for i=1:t
    eval(['n(i)=length(locsy' num2str(i) ');']);
    nmin=min(n);
end;

%% put locsy# variable into a matrix and make the periods into a matrix
for i=1:t
    eval(['locsmx(1:nmin,i)=locsy' num2str(i) '(1:nmin);']);
    for j=1:nmin-1
        taumx(j,i)=2*(locsmx(j+1,i)-locsmx(j,i));
    end;
end;
%%%% create time variable for plotting
time=2*locsmx(1:nmin-1,:);

figure
plot(time,taumx)
%% Calculate phase shift phi
expdphi=zeros(nmin-1,t-1);
for i=1:t-1
    for j=1:nmin-1
        % %         if i<t
        % %             %         phimx(j,i)=min([abs(locsmx(j,i+1)-locsmx(j,i))/taumx(j,i), abs(1-abs(locsmx(j,i+1)-locsmx(j,i))/taumx(j,i)), abs(2-abs(locsmx(j,i+1)-locsmx(j,i))/taumx(j,i))]);
        % %             %         else phimx(j,i)=min([abs(locsmx(j,i)-locsmx(j,1))/taumx(j,i), abs(1-abs(locsmx(j,i)-locsmx(j,1))/taumx(j,i)), abs(2-abs(locsmx(j,i)-locsmx(j,1))/taumx(j,i))]);
        % %             phimx(j,i)=min(abs(locsmx(j,i+1)-locsmx(j,i))/taumx(j,i), abs(2-abs(locsmx(j,i+1)-locsmx(j,i))/taumx(j,i)));
        % %         else phimx(j,i)=min(abs(locsmx(j,i)-locsmx(j,1))/taumx(j,i), abs(2-abs(locsmx(j,i)-locsmx(j,1))/taumx(j,i)));
        % % %             phimx(j,i)=(locsmx(j,i+1)-locsmx(j,i))/taumx(j,i);%raw data
        % % %         else phimx(j,i)=(locsmx(j,i)-locsmx(j,1))/taumx(j,i);
        % %         end;
        %         dphi(j,i)=min(abs(locsmx(j,i)-locsmx(j,1))/(taumx(j,i)/2),abs(2-abs(locsmx(j,i)-locsmx(j,1))/(taumx(j,i)/2)));
        %         dphi(j,i)=abs(locsmx(j,i)-locsmx(j,1))/(taumx(j,i)/2);
%%%%%%%%%to avoid ambiguity when rebuilding stplot from dphi
        dphi=round(locsmx(j,i+1)-locsmx(j,1))/(taumx(j,i+1)/2);
        while dphi<0 || dphi>1
            if dphi<0
                dphi=dphi+1;
            else
                dphi=dphi-1;
            end;
        end;
        expdphi(j,i)=dphi;
    end;
end;
figure
plot(time(:,2:t),expdphi)
