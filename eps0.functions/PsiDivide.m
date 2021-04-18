function [psiw,psib]=PsiDivide(psi,lw,zArr) 

%Store the component of the wave-function inside the well and barrier
%respectively. 

for ia=1:size(psi,2)
temp=psi(:,ia);

psiw(:,ia)=(temp(abs(zArr)<=lw/2));
psib(:,ia)=(temp(abs(zArr)>lw/2));
end

end