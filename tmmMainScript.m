clear,clc
addpath('tmmScripts');

%Specify polarization (either "tm" or "te") and the array of incident
%angles and the height of the structure. In "Ballistic Metamaterials" the
%total height of the periodic structure was 4um. This code recreates the
%absorption from HMM3 and uses the 4um total height of the stack. 

polar='tm'; %polarization of incident light, either "tm" or "te".
anglist=(5:1:75); %angle of incident light. 
strucHt=4; %height of composite measured in microns. 

%Load in the permittivity data (from running epsMainScript.m). Remember,
%run epsMainScript.m FIRST!

epsMat=load('epsMat.lay=9.5nm.mat'); 

%Load in parameters from running "epsMainScript.m" FIRST.
lamlist=epsMat.lamPlot;
epsXYArr=epsMat.epsParl;
epsZZArr=epsMat.epsPerp;


%Obtain optical coefficients.
[rTot,tTot,aTot]=spectralProps(lamlist,anglist,strucHt,epsXYArr,epsZZArr,polar); 

%rTot - reflection from composite. 
%tTot - transmision through composite. 
%aTot - absorption of composite. 

%% plotting 
figure(1)
surf(lamlist,anglist,aTot,'edgecolor','none'); view(2);
colormap hot; xlim([2 14]); ylim([5 75]); 
xlabel('Wavelength [um]'); ylabel('Inc. Angle [deg.]'); 
colorbar; title(colorbar,'A'); caxis([0 1]); title('HMM3');
set(gca,'fontsize',18)

clear cT