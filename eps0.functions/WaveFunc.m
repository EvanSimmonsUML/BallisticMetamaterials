function psi=WaveFunc(bndst,kw,kb,zArr,lw,data)

%(9/28/2018) : This function calculates the wave-function of a finite
%SYMMETRIC square well. NOTE: lw must be in nm for this to work. 

lw=lw*data.nm;
zArr=zArr*data.nm;

for ia=1:length(bndst)
    if mod(ia,2)~=0
        alpha=kw(ia)*lw/2;
        beta=kb(ia)*lw/2;
        
        a=cos(alpha)/(exp(-beta));
        
        for ib=1:length(zArr)
            if zArr(ib)<=-lw/2
                psi(ib,ia)=a*exp(kb(ia)*zArr(ib));
            elseif zArr(ib)>-lw/2 && zArr(ib)<lw/2
                psi(ib,ia)=cos(kw(ia)*zArr(ib));
            else
                psi(ib,ia)=a*exp(-kb(ia)*zArr(ib));
            end
        end
        
    else
        
        alpha=kw(ia)*lw/2;
        beta=kb(ia)*lw/2;
        
        b=sin(alpha)./(exp(-beta));
        
        for ic=1:length(zArr)
            if zArr(ic)<=-lw/2
                psi(ic,ia)=-b*exp(kb(ia)*zArr(ic));
            elseif zArr(ic)>-lw/2 && zArr(ic)<lw/2
                psi(ic,ia)=sin(kw(ia)*zArr(ic));
            else
                psi(ic,ia)=b*exp(-kb(ia)*zArr(ic));
            end
        end
        
    end
    
end