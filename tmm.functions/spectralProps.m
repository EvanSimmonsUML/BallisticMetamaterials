function [rTot,tTot,aTot]=spectralProps(lamlist,anglist,strucHt,epsXYArr,epsZZArr,polar)

%planar TMM script used to calculate reflection,transmission, and
%absorption using optically thick substrate analysis. 

c=3e8;
rTot=zeros(length(anglist),length(lamlist)); 
tTot=rTot; 

ainc=1;

for iang=1:length(anglist)
    for ilam=1:length(lamlist)
        zFBArr=[0 strucHt];
        
        epsSub=subEps(lamlist(ilam));

        omg=2*pi*c/(lamlist(ilam));
        
        %air-hmm-sub
        epsXYTot=[1 epsXYArr(ilam) epsSub];
        epsZZTot=[1 epsZZArr(ilam) epsSub];
        
        kxArr=omg*sqrt(epsXYTot(1,1))*sind(anglist(iang))/c;
        [rF,tF]=tmmFun(omg,kxArr,epsXYTot,epsZZTot,zFBArr,ainc,polar);
        
        %sub-hmm-air
        epsXYTot=[epsSub epsXYArr(ilam) 1];
        epsZZTot=[epsSub epsZZArr(ilam) 1];
        
        angRef=asind(sind(anglist(iang))/sqrt(epsSub));
        
        kxArr=omg*sqrt(epsXYTot(1,1))*sind(angRef)/c;
        [rB,tB]=tmmFun(omg,kxArr,epsXYTot,epsZZTot,zFBArr,ainc,polar);
        
        %air-sub
        epsXYTot=[1 epsSub];
        epsZZTot=epsXYTot;
        zFBArr=0;
        
        kxArr=omg*sqrt(epsXYTot(1,1))*sind(anglist(iang))/c;
        [rS,tS]=tmmFun(omg,kxArr,epsXYTot,epsZZTot,zFBArr,ainc,polar);
        
        rTot(iang,ilam)=rF+(tF*tB*rS)/...
            (1-rS*rB);
        
        tTot(iang,ilam)=(tF*tS)/...
            (1-rS*rB);
    end
end

aTot=1-rTot-tTot;
end