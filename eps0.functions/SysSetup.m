%% SysSetup %%

%Initialize all SI conversions and material parameters. Accepts user
%defined material parameters and returns structure of all material data and
%constants. 

function [data]=SysSetup(v0Ht,m0W,m0B,npParam,epsW,epsB)

m_e=9.11e-31;
h_bar=6.63e-34/2/pi;
e=1.602e-19; 
eps0=8.85e-12;
c=3e8;

ev2joule=e;
joule2ev=1/e;
cm2m=1e2;

nm=1e-9;
um=1e-6;

%Structure of Constants% 
data.m_e=m_e;
data.h_bar=h_bar;
data.e=e;
data.eps0=eps0;
data.ev2joule=ev2joule;
data.joule2ev=joule2ev;
data.cm2m=cm2m;
data.nm=nm;
data.um=um;
data.c=c;
data.mW=m0W*data.m_e;
data.alpha=npParam;
data.v0=v0Ht;
data.mB=m0B*data.m_e;
data.epsW=epsW;
data.epsB=epsB;
