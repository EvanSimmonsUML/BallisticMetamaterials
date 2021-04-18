function [epsMet]=epsCalcScript(layW,layB,kPMax,w,dmpW,v0Ht,data)

%This script returns the permittivity of the highly doped layer in the
%system using a finite (isolated) quantum well model. 
zpnts=2e3; 
zArr=linspace(-(layW+layB)/2,(layW+layB)/2,zpnts);
dz=(zArr(end)-zArr(1))/(length(zArr)-1);

kPArr=linspace(0,kPMax,1e3)*pi/layW/data.nm;
dkP=(kPArr(end)-kPArr(1))/(length(kPArr)-1);

aeta=0;

for ikp=1:length(kPArr)
    tic=0;
    
    if ikp==1
        %Calculate Bound Energies.
        [bndst,mw,mb]=EigenE(layW,data);
        
        [kb,kw]=kvec(mw,mb,bndst,data);
        
        %Calculate Wave Functions
        psi=WaveFunc(bndst,kw,kb,zArr,layW,data);
        
        %Normalize Wave Functions
        psi=WaveNorm(psi,dz,data);
        
        %Calculate Fermi Energy
        fermi=fermiScript(bndst,data,layW);
    end
    
    
    %Find portion of Wave Function in each Material
    [psiw,psib]=PsiDivide(psi,layW,zArr);
    [Pw,Pb]=ProbDensity(psiw,psib,bndst,dz,data);
    
    mw0=effMScript(bndst,data);
    mb0=data.mB;
    
    %In-Plane Effective Mass
    mxy=Pw.*mw0+Pb.*mb0;
    
    %Calculate subband dispersion
    for ibnd=1:length(bndst)
        Esub(ibnd,ikp)=bndst(ibnd)+...
            data.h_bar^2*kPArr(ikp)^2/(2*mxy(ibnd))*data.joule2ev;
        
        if fermi>v0Ht
            disp(['E_f = ',num2str(fermi),'. Exceeds barrier height!']);
            break;
        end
        
        if Esub(ibnd,ikp)<=fermi && Esub(ibnd,ikp)>0
            tic=tic+1;
        end
        
        
    end
    
    %Calculate Oscillator Strength and Transition Frequencies
    [~,fd{ikp}]=Freq(Esub(:,ikp),bndst,psi,zArr,dz,data);
    [wt{ikp},~]=FreqAdjust(Esub(:,ikp),fermi,psi,data);
    aeta=aeta+tic*kPArr(ikp)*dkP;
    
    if Esub(1,ikp)>fermi
        break;
    end
end

aeta=data.n_dop/aeta;

Esub(Esub==0)=nan;

D=@(wtx,dm) wtx^2-w.^2-dm*wtx*1i*w;
chi=zeros(1,length(w));

%Calculate permittivity
for ih=1:ikp
    wtMat=wt{ih};
    fdMat=fd{ih};
    
    wp=sqrt(aeta*kPArr(ih)*dkP*data.e^2./effMScript(fermi,data)/data.eps0/data.epsW);
    
    for ii=1:size(wtMat,1)
        for ij=1:size(wtMat,2)
            if wtMat(ii,ij)>0
                chi=chi+wp.^2*fdMat(ii,ij)./D(wtMat(ii,ij),dmpW);
            end
        end
    end
end
epsMet=data.epsW*(1+chi);
end