shot = 78071;
t = 1.0;
[L,~,LY] = liu('tcv',shot,t);

%%
Lfge = fge('tcv',LY.shot); % get FGE (another equilibrium solver)
LXfge = meqxconvert(L,LY,Lfge); % convert LX data structure based on L into LX2 structure based on L2 
LXfge.bpD = LXfge.bp; % Beta poloidal (per domain)

%%
for i = 1 : numel(LY.t)
  %  FGE linearization around initial equilibrium guess LX
  Lfge_(i) = fgel(Lfge,meqxk(LXfge,i)); 
  % fgeeig Calculate unstable eigenvalues of system using fgess
  ee(i)    = fgeeig(Lfge_(i)); 
end