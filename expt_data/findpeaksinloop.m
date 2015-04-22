%created Sep 03, 2009 by Ning Li for running findpeaks in a loop of files

%% locate peaks
t = 3; %total number of profiles (droplets)
name = 'sa'; %profile file name
dis = 60; %minpeakdistance

for i=1:t
    
    s=[ name num2str(i) '.txt'];
    load (s);
    eval(['x' num2str(i) '=' name num2str(i) '(:,1);']);
    eval(['y' num2str(i) '=' name num2str(i) '(:,2);']);
    eval(['[pksy' num2str(i) ',locsy' num2str(i) ']= findpeaks ( y' num2str(i) ', ''minpeakdistance'',  dis );']);
    
    figure(i)
    eval(['plot(x' num2str(i) ',y' num2str(i) ',locsy' num2str(i) ',pksy' num2str(i) ', ''marker'' , ''o'' )']);
    
end

%% calculate periods
for i=1:t
    eval(['[m,n]=size(locsy' num2str(i) ');']);
    for j=1:n-1
        eval(['tau_' num2str(i) '(j) = locsy' num2str(i) '(j+1)-locsy' num2str(i) '(j);']);
    end
end