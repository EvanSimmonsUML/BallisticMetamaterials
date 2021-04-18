%% FreqAdjust %%

%Calculate the actual transitions using fermi energy%

function [wt,EdifList]=FreqAdjust(Esub,fermi,psi,data)

wt=zeros(length(Esub),length(Esub));
EdifList=wt;

for ia=1:size(psi,2)
    for ib=1:size(psi,2)
        %       Generalized Transition Frequency Matrix%
        %       ROWS: Initial state
        %       COLUMNS: Final State
        %       omgfreq(x,y)->x is initial, y is final
        
        if Esub(ib)>fermi && Esub(ia)<=fermi && (Esub(ib))~=0 && (Esub(ia))~=0 && Esub(ia+1)>fermi 
            wt(ia,ib)=(Esub(ib)-Esub(ia))*data.ev2joule/data.h_bar;
            EdifList(ia,ib)=(Esub(ib)-Esub(ia));
        end
    end
end

end
