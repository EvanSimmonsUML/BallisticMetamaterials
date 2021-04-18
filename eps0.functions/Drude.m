%% Drude %%
%Calculate Classical Dielectric Function. 
function epsD=Drude(w,wp,T,epsInf)
epsD=epsInf*(1-wp^2./(w.^2+1i*w./T));
end