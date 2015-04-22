%created Sep 03, 2009 by Ning Li for running findpeaks in a loop of files

close all
clear
%% find peaks
t = 1; %total number of profiles (droplets)
name = 'c'; %profile file name
dis = 100; %minpeakdistance

for i=1:t
    
    s=[ name num2str(i) '.txt'];
    load (s);
    eval(['x' num2str(i) '=' name num2str(i) '(:,1);']);
    eval(['y' num2str(i) '=' name num2str(i) '(:,2);']);
    
    eval(['frame=length(x' num2str(i) ');']);
    %eval(['xmx(1:frame,i)= x' num2str(i) '(1:frame);']);
    eval(['ymx(1:frame,i)= y' num2str(i) '(1:frame);']);
    
    
    % deal with flat peak (same maxima for two neigboring points)
    for j=1:frame-1
        if ymx(j,i) == ymx(j+1,i)
            ymx(j,i)= ymx(j,i)-0.1;
        end
    end
    
    %eval(['xn' num2str(i) '= xmx(1:frame,i);']);
    eval(['yn' num2str(i) '(1:frame,1) = ymx(1:frame,i);' ]);
    %now it should be fine
    
    eval(['[pksy' num2str(i) ',locsy' num2str(i) ']= findpeaks ( yn' num2str(i) ', ''minpeakdistance'',  dis );']);
    
    figure(i)
    eval(['plot(x' num2str(i) ',yn' num2str(i) ',locsy' num2str(i) '-1,pksy' num2str(i) ', ''marker'' , ''o'' )']);
    
end

%% calculate periods
figure
hold all
for i=1:t
    eval(['[m,n]=size(locsy' num2str(i) ');']);
    eval(['tau_' num2str(i) '=zeros(1,m-1);']);

    for j=1:m-1
        eval(['tau_' num2str(i) '(j) = 4*(locsy' num2str(i) '(j+1)-locsy' num2str(i) '(j));']);
    end
    eval(['plot(tau_' num2str(i) ')']);
end
hold off

tauavg=mean(tau_1)