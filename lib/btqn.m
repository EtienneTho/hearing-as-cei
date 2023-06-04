function info = btqn(wm, criterion)
% function info = btqn(wm)
% 
% btqn uses a quasi-Newton-accelerated MM algorithm to fit the
% Bradley-Terry model.
%
% The input wm is an nxn matrix such that wm(i,j) is the number of times 
% team i beat team j.
%
% The output is an nx1 vector of parameter estimates of "team skill".
%
% This algorithm does not contain any checks for lack of convergence;
% it is assumed that for any partition of the teams into sets A and B,
% at least one team from A beats at least one team from B at least once.

n=size(wm,1);
% flops(0);

% The parameter vector gamma is constrained to keep its last
% entry equal to 0; thus, nmo stands for "n minus one".
nmo=n-1;
gamma=zeros(nmo,1); % initial value:  All teams equal
iteration=0;
change=realmax;
gm=wm(1:nmo,:)+wm(:,1:nmo)';
wins=sum(wm(1:nmo,:),2);
gind=(gm>0);
y=zeros(nmo,n);
z=y;
M=zeros(nmo,nmo);
newchange=0;

while norm(change) > criterion  % Simplistic convergence criterion 
	iteration = iteration + 1;
	pi=exp(gamma);
	pi(n)=1;
		pius=pi(:,ones(n,1));
		piust=pius(:,1:nmo)';
		pius=pius(1:nmo,:);
		y(gind) = gm(gind) ./ (pius(gind)+piust(gind));
		z(gind)=pius(gind).*y(gind);
		dL=wins - sum(z,2);
		newpi=wins ./ sum(y,2);
	newgamma=log(newpi);
	oldchange=newchange;
	newchange=newgamma-gamma;
	if iteration>1
		g=olddL-dL;
		s=-change+oldchange-newchange;
		s_plus_Mg=s+M*g;
		M=M-s_plus_Mg * s_plus_Mg' / (g' * s_plus_Mg);
	end	
	olddL=dL;
	change=newchange+M*dL;
    if isnan(change)
        break
    end
	gamma=gamma+change;		
end
% Iterations = iteration
% Floating_point_operations=flops
info=exp(gamma);
info(n)=1;

