%{
Visualize Hilbert transform for sample fMRI timeseries for the methods
description.
Adapts code from: 
https://github.com/juanitacabral/LEiDA_Psilocybin
Output:
fMRI_hilbert_1.jpg Hilbert transform plot for one ROI.
fMRI_hilbert_2.jpg Hilbert transform plot for a second ROI. 
%}

%Define command line arguments.
function [] = r_hilbert_plot(subject)
disp(append('Doing: ',subject))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Load BOLD from one subject and one region.
TR = 0.72;
inpath = append('../outputs/r_meants/',subject,'/');
infile = append(inpath,'demean_postproc_rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');
tc = readmatrix(infile);
croi = 2;
signal = tc(croi,:);

%Hilbert transform.
Complex_BOLD = hilbert(signal);
Theta=angle(Complex_BOLD(1:40));

%Plot.
global ColorOrder, ColorOrder='kkk';
figure
plot3(0:TR:(length(Theta)-1)*TR,sin(Theta),cos(Theta),'k','LineWidth',1,'Color','red');
xlabel('Time (seconds)')
ylabel('Imag')
zlabel('Real')
grid on
ylim([-3 3])
zlim([-3 3])
hold on
P1=[0:TR:(length(Theta)-1)*TR;zeros(size(Theta));zeros(size(Theta))]';
P2=[0:TR:(length(Theta)-1)*TR;sin(Theta);cos(Theta)]';
daspect([3.5 1 1])
arrow3(P1,P2,'o',.5)
arrow3([-10 0 0],[length(Theta)*TR+10 0 0],'o')
xlim([-10 30])
hold on
plot3(0:TR:(length(Theta)-1)*TR,3*ones(size(Theta)),cos(Theta),'r');
plot3(0:TR:(length(Theta)-1)*TR,sin(Theta),-3*ones(size(Theta)),'r:');
P1(:,1)=-10;
P2(:,1)=-10;
arrow3(P1,P2,'o',0)

%Orient.
[caz,cel] = view;
view(-(caz-20),cel)

%Save.
outfile = '../outputs/outcollect/fMRI_hilbert_1.jpg';
exportgraphics(gcf,outfile,'Resolution',720);

%Load BOLD from one subject and another region.
croi = 1;
signal = tc(croi,:);

%Hilbert transform.
Complex_BOLD = hilbert(signal);
Theta=angle(Complex_BOLD(1:40));

%Plot.
global ColorOrder, ColorOrder='kkk';
figure
plot3(0:TR:(length(Theta)-1)*TR,sin(Theta),cos(Theta),'k','LineWidth',1,'Color','blue');
xlabel('Time (seconds)')
ylabel('Imag')
zlabel('Real')
grid on
ylim([-3 3])
zlim([-3 3])
hold on
P1=[0:TR:(length(Theta)-1)*TR;zeros(size(Theta));zeros(size(Theta))]';
P2=[0:TR:(length(Theta)-1)*TR;sin(Theta);cos(Theta)]';
daspect([3.5 1 1])
arrow3(P1,P2,'o',.5)
arrow3([-10 0 0],[length(Theta)*TR+10 0 0],'o')
xlim([-10 30])
hold on
plot3(0:TR:(length(Theta)-1)*TR,3*ones(size(Theta)),cos(Theta),'k:');
plot3(0:TR:(length(Theta)-1)*TR,sin(Theta),-3*ones(size(Theta)),'k:');
P1(:,1)=-10;
P2(:,1)=-10;
arrow3(P1,P2,'o',0)

%Orient.
[caz,cel] = view;
view(-(caz-20),cel)

%Save.
outfile = '../outputs/outcollect/fMRI_hilbert_2.jpg';
exportgraphics(gcf,outfile,'Resolution',720);
end
