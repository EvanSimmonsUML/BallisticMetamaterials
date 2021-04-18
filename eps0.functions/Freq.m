%% Freq %%

%Calculate the transition frequency and transition dipole.
function [wt,fd]=Freq(Esub,bndst,psi,zArr,dz,data)

wt=zeros(length(Esub),length(Esub));
fd=wt;

for ia=1:size(psi,2)
    for ib=1:size(psi,2)
        %       Generalized Transition Frequency Matrix%
        %       ROWS: Initial state
        %       COLUMNS: Final State
        %       omgfreq(x,y)->x is initial, y is final  
            wt(ia,ib)=(Esub(ib)-Esub(ia))*data.ev2joule/data.h_bar;
            %Calculate all possible transition energies. 
            if wt(ia,ib)>0
            dpole=Dipole(psi(:,ia),psi(:,ib),zArr,dz,data);
            
            %NP fd - Av. Energy%
            fd(ia,ib)=2*effMScript(mean([bndst(ia),bndst(ib)]),data)*wt(ia,ib)*abs(dpole)^2/data.h_bar;
            end
    end
end
end