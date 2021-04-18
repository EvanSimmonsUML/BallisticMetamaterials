function [rC,tC,pSzArr,nSzArr,...
    pEflux,nEflux] = tmmFun(omg,kxArr,eps1j,eps2j,zArr,ainc,p)

%explicit TMM calculation for anisotropic planar media.
mu0=4*pi*1e-7; c=3e8;

kzArr=zeros(1,length(eps1j));
mMat=cell(1,length(eps1j));
tMat=cell(1,length(eps1j));
apArr=zeros(1,length(eps1j));
anArr=apArr;
pSzArr=apArr;
nSzArr=apArr;

if strcmp(p,'tm')==true%TM
for ii = 1:length(eps1j)
    kzArr(ii)=sqrt((omg/c)^2*eps1j(ii)-((eps1j(ii)/eps2j(ii))*kxArr^2));
    mMat{ii} = [kzArr(ii)*eps2j(ii), -kzArr(ii)*eps2j(ii); omg*eps1j(ii)*eps2j(ii)/c/c,omg*eps1j(ii)*eps2j(ii)/c/c];
end 

else
    for ii = 1:length(eps1j)
    kzArr(ii) = sqrt(eps1j(ii)*omg^2/c^2 - kxArr^2);
    mMat{ii} = [1 1;-kzArr(ii)/omg/mu0 kzArr(ii)/omg/mu0];
    end
end 
%---------------------------------------------------%
       %Calculations of PSI and T matrices%
%---------------------------------------------------%

TMtot = eye(2);
for jj = 1:length(eps1j)-1
    PSI = diag(exp([1i*kzArr(jj)*zArr(jj);-1i*kzArr(jj)*zArr(jj)]));
    PSI_1 = diag(exp([1i*kzArr(jj+1)*zArr(jj);-1i*kzArr(jj+1)*zArr(jj)]));
    
    tMat{jj} = (mMat{jj+1}*PSI_1)\(mMat{jj}*PSI);
    
    TMtot = tMat{jj}*TMtot;
    
end 
%---------------------------------------------------%
           %Calculations of amplitudes%
%---------------------------------------------------%

aref = -TMtot(2,2)\TMtot(2,1)*ainc;

apArr(1) = ainc;
anArr(1) = aref;

aj = [ainc;aref];

for il = 1:length(eps1j)-1
    aj1 = tMat{il}*aj;
    apArr(il+1) = aj1(1);
    anArr(il+1) = aj1(2);
    aj = aj1;
end 

%---------------------------------------------------%
         %Poynting Vector Calculations%
%---------------------------------------------------%
if strcmp(p,'tm')==true
for jl = 1:length(eps1j)
    pEfield = [kzArr(jl)*eps2j(jl);0;-eps1j(jl)*kxArr];
    nEfield = [-kzArr(jl)*eps2j(jl);0;-eps1j(jl)*kxArr];
    pHfield = [0;omg*eps1j(jl)*eps2j(jl)/c/c;0];
    nHfield = pHfield;
    
    pSzArr(jl) = .5*real(dot([0;0;1],cross(pEfield,conj(pHfield))));
    nSzArr(jl) = .5*real(dot([0;0;1],cross(nEfield,conj(nHfield))));
    
end 

else
    for jl = 1:length(eps1j)
    pEfield = [0;1;0];
    nEfield = pEfield;
    pHfield = [-kzArr(jl);0;kxArr]/omg/mu0;
    nHfield = [kzArr(jl);0;kxArr]/omg/mu0;
    
    pSzArr(jl) = .5*real(dot([0;0;1],cross(pEfield,conj(nHfield))));
    nSzArr(jl) = .5*real(dot([0;0;1],cross(nEfield,conj(pHfield))));
    end
end 

pEflux = pSzArr.*abs(apArr).^2;
nEflux = nSzArr.*abs(anArr).^2;

rC=-nEflux(:,1)./pEflux(:,1);
tC=pEflux(length(eps1j))./pEflux(:,1);
end