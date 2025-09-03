shot = 78071;
t = 1.0;
[L,~,LY] = liu('tcv',shot,t);

%%
Lfge = fge('tcv',LY.shot);
LXfge = meqxconvert(L,LY,Lfge);
LXfge.bpD = LXfge.bp;

%%
for i = 1 : numel(LY.t)
  Lfge_(i) = fgel(Lfge,meqxk(LXfge,i));
  ee(i)    = fgeeig(Lfge_(i));
end