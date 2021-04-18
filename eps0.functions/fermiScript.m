function Ef=fermiScript(bndst,data,lw)

%Calculate the fermi energy using the bound states. Ensure convergence of
%the energy. 

%Conversions%
lw=lw.*data.nm;
%Constants%
zeta=data.n_dop*lw*pi*data.h_bar^2/data.mW;
zeta=zeta*data.joule2ev;

%Fermi Calculations%
if length(bndst)==1
    A=data.alpha;
    B=1;
    C=-1*(zeta+(data.alpha)*(bndst(1).^2)+(bndst(1)));
    Ef=(-B+sqrt(B^2-4*A*C))/2/A;
else
    for sN=1:length(bndst)
        A=data.alpha;
        B=1;
        C=-(1/sN)*(zeta+(data.alpha)*sum(bndst(1:1:sN).^2)+sum(bndst(1:1:sN)));
        Ef=(-B+sqrt(B^2-4*A*C))/2/A;
        
        if   sN~=length(bndst)
            if Ef<bndst(sN+1)
                break;
            elseif sN==length(bndst)
                break;
            end
        end
    end
end
end