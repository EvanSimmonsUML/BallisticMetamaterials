%% kvec %%

%Calculate well & barrier k vector.

function [kb,kw]=kvec(mw,mb,bndst,data)

kb=sqrt((2.*mb/data.h_bar./data.h_bar).*(data.v0-bndst).*data.ev2joule);
kw=sqrt(2.*bndst.*data.ev2joule.*mw./data.h_bar./data.h_bar);

end