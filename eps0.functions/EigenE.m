function [bndst,effmW,effmB]=EigenE(lw,data)
%Calculate bound state energies using nonparabolic effective mass approach.
%
tol=1;
eArr=linspace(0,data.v0,1e6); %eV

lw=lw*data.nm;
vOff=data.v0*data.ev2joule;
eArrjl=eArr.*data.ev2joule;


for ia=1:length(eArr)
    
    kw=(1/data.h_bar).*sqrt(2.*effMScript(eArr(ia),data).*eArrjl(ia));
    kb=(1/data.h_bar).*sqrt(2*data.mB.*(vOff-eArrjl(ia)));
    
    alpha=kw.*lw./2;
    
    zeroFun(ia)=(tan(alpha)-kb.*effMScript(eArr(ia),data)./(kw.*...
        data.mB)).*(cot(alpha)+(kb.*effMScript(eArr(ia),data))./...
        (kw.*data.mB));
    
    %Boundary Condition: http://aip.scitation.org/doi/pdf/10.1063/1.3457787 
    kw=0;
    kb=0;
end


for ib=1:length(zeroFun)-1
    holder=zeroFun(ib)*zeroFun(ib+1);
    if holder<0 && zeroFun(ib)<tol
        bndst(ib)=(eArr(ib+1)+eArr(ib))/2;
    end
end

bndst=bndst(bndst~=0);
effmW=effMScript(bndst,data);
effmB=data.mB;
end
