%Calculate dipole integral.
function [dpole]=Dipole(Ipsi,Fpsi,zArr,dz,data)

zArr=zArr*data.nm;
dz=dz*data.nm;
dpole=0;

dpole=sum(conj(Fpsi).*zArr'.*Ipsi).*dz;
end
