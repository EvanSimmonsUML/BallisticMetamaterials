%% mEff %%
%Calculate effective mass of InGaAs at E.
function effm=effMScript(E,data)
effm=data.mW*(1+data.alpha*E);
end