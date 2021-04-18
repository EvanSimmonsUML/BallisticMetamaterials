function [pw,pb]=ProbDensity(psiw,psib,bndst,dz,data)

%Probability (quantum mechanical using wave functions) of electron being
%inside well and barrier respectively. 

for ia=1:size(psiw,2)
    temp1j=psiw(:,ia);
    temp2j=psib(:,ia);
    
    pw(ia)=(temp1j)'*temp1j.*dz.*data.nm;
    pb(ia)=(temp2j)'*temp2j.*dz.*data.nm;
end
end