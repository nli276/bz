% S coupling function contour for phase diagram
close all
clear
%%
a=transpose(100e-6:20e-6:180e-6);
% ma=10e-3:10e-3:100e-3;
b=(100e-6:20e-6:180e-6);
% a=transpose(5e-6:5e-6:25e-6);%m
% b=(5e-6:5e-6:25e-6);
% a=transpose(2e-6:2e-6:50e-6);
% b=(2e-6:2e-6:50e-6);
% a=transpose(1e-6:1e-6:10e-6);%m
% b=(1e-6:1e-6:10e-6);
ma=0.4;%M

pB=2.5;
D=1e-9;%m^2/s
% S1=pB*D./(2*(a.^2)*70*ma);

S2=pB*D/(ma*70)./(a*b);%symmetric
S1=pB*D/(ma*70)./(a*b+repmat(a.^2,1,length(b)));%Marcin, asymmetric
% S1=pB*D/(ma*70)./(a*b+repmat(pB*a.^2,1,length(b)));%New Vlad

% S0=pB*D/(10+29*ma)./(a*b);
% S0=pB*D./((a.^2)*(10+29*ma));%turns out S1 fits phase diagram better
%%
hold all
axis=1:10;
% axis=1:5;
% axis=1:20;

[C,h]=contour(axis,axis,S1,[0.05, 0.1, 0.2]);%'LevelStep',0.01);%
% axis([xmin xmax ymin ymax])
% axis([0.2 0.4 1.5e-4 2.5e-4])
clabel(C,'manual');
xlabel('oil size (\mum)')
% xlabel('[MA] (mM)')
ylabel('drop size (\mum)')
% title('S')