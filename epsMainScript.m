clear,clc
addpath('eps0Scripts');

%% Section 1: Unit Cell Illustration

%Script calculates the dielectric tensor components for a highly
%anisotropic composite using the methodology outlined in the paper
%"Ballistic Metamaterials" (https://doi.org/10.1364/OPTICA.402891). The
%unit cell is as follows, 
%
%
%                     v0Ht          v0Ht
%       ---------------|             |---------------  
%                      |             |                
%                      |             |               
%          Barrier     |    Well     |    Barrier               
%          m0B,layB    |(Doped Layer)|    m0B,layB           
%                      |             |               
%          (AlInAs)    |  (InGaAs)   |    (AlInAs)            
%                      |   npParam   |               
%                      |  layW,m0W   |
%                      |-------------|
%
% Figure shows relevant variables (used in script) and setup for isolated
% finite quantum well. The figure illustrates the barrier and well material and 
% the energy landscape of the system. Note: codes do not incorporate superlattice effects
% (Kronig-Penney etc.) as discussed in Manuscript.
%

%% Section 2: Material Parameters%%
%System defined for HMM3: Heavily doped InGaAs (layer thickness of 9.5nm)
%separated by undoped AlInAs (layer thickness of 9.5nm too) (schematic
%above). 

%Well Material: InGaAs ; Barrier Material: AlInAs

layW=9.5; %thickness of the well material in nanometers.
layB=layW; %thickness of the barrier material in nanometers.
dopNum=1.78e+19; %electron concentration inside well material in cm^-3.
v0Ht=.52; %Barrier height from the heterojunction formed in eV. 
epsW=12.15; %Background permittivity of the well material.
epsB=10.23; %Background permittivity of the barrier material.

kPMax=25; %Max lateral wavenumber for in-plane energy calculations.
          %Converted into units of inverse meters (standard for energy dispersion
          %calculations).  

m0W=0.0427; %Effective Mass Coefficient at the base of the conduction band for the well
            %material. Actual Effective mass is m0W*m_e with m_e=9.11e-31 kg. 

m0B=0.075; %Effective Mass Coefficient at the base of the conduction band
           %for the barrier material. Actual effective mass in the barrier is
           %m0B*m_e with m_e=9.11e-31kg.
           
           
npParam=1.24; %non-parabolicity parameter for the well material. Given in units of 1/eV this
              %value is often given as the inverse of the band-gap energy
              %of the material. 


data=SysSetup(v0Ht,m0W,m0B,npParam,epsW,epsB); %Setup list of various constants in SI units. 

data.n_dop=dopNum.*data.cm2m^3; %Convert doping to SI units. 

%% Section 3: Optical Parameters%%
%Drude Parameters for in-plane permittivity component.

wpBlk=8.5; %Plasma wavelength associated with plasma frequency of material. Given in
           %microns. 
dmpBlk=1e13; %scattering rate for bulk material for Drude Model. Given in Hertz. 

lamMin=2; %Minimum wavelength for permittivity calculations - given in microns. 
lamMax=14; %Maximum wavelength for permittivity calculations - given in microns. 

dmpW=0.1; %Transition strength coefficient. Commonly set to 10% of the resonance
          %frequency of an optical transition. 

%% End User Input: Begin code setup. 
%Wavelength range and frequency range. 
lamlist=linspace(lamMin,lamMax,1200)*data.um; lamPlot=lamlist*1e6;
w=2*pi*data.c./lamlist;

%Calculate out-of-plane permittivity of highly doped layer. 
epsMet=epsCalcScript(layW,layB,kPMax,w,dmpW,v0Ht,data);

%Calculate in-plane permittivity of highly doped layer.
epsBlk=Drude(w,2*pi*data.c./(wpBlk*data.um),1/dmpBlk,data.epsW);

%Calculate permittivity components of total composite.
epsParl=(layW*epsBlk+layB*data.epsB)./(layW+layB);
epsPerp=((layW+layB).*epsMet.*data.epsB)./(layW*data.epsB+...
    layB.*epsMet);

%% plotting
figure(1)
plot(lamPlot,real(epsParl),lamPlot,real(epsPerp),'linewidth',2);
xlabel('Wavelength [um]'); ylabel('\epsilon_{EMT}'); set(gca,'fontsize',18);
xlim([lamMin lamMax]); title('HMM3'); 
grid on; box on; 
legend('\Re (\epsilon_{||})','\Re (\epsilon_\perp)','location','northwest');

%Save data for spectral calculations.
save(['epsMat.lay=',num2str(layW),'nm.mat'],'lamPlot','epsParl','epsPerp');