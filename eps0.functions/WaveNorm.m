%% PsiNorm %%

%Calculate Normalized Wave Function.

function [PsiNorm,A]=WaveNorm(PsiArr,dz,data)
dz=dz*data.nm;

for ia=1:size(PsiArr,2)
    PsiPsi=PsiArr(:,ia)'*PsiArr(:,ia);
    PsiInt(ia)=PsiPsi*dz;
    A(ia)=1/sqrt(PsiInt(ia));
    PsiNorm(:,ia)=A(ia).*PsiArr(:,ia);
end