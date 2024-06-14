function [k,p] = dlqri(a,b,q,r,del)
% DLQRI  Iterative linear quadratic regulator design.
%
% [k,p] = dlqri(a,b,q,r,del)
%
% a,b,c,d: State space form.
% q,r: State and input weights.
% del: Convergence tolerence (default=1e-8).
%
% k: Optimal feedback gain matrix.
% p: Final P matrix.
%
% Based on Control Toolbox DLQR, but modified into an iterative
% form to deal with singular A matrices.

% James Taylor, 09/08/1999
% Wlodek Tych, 02/02/1990

error(nargchk(4,5,nargin));
error(abcdchk(a,b));
if nargin<5, del=1e-8; end
[m,n] = size(a);
[mb,nb] = size(b);
[mq,nq] = size(q);
if (m ~= mq) | (n ~= nq) 
	error('A and Q must be the same size')
end
[mr,nr] = size(r);
if (mr ~= nr) | (nb ~= mr)
	error('B and R must be consistent')
end

% check if q is positive semi-definite and symmetric
if any(eig(q) < 0) | (norm(q'-q,1)/norm(q,1) > eps)
	error('Q must be symmetric and positive semi-definite')
end
% check if r is positive definite and symmetric
if any(eig(r) <= 0) | (norm(r'-r,1)/norm(r,1) > eps)
	error('R must be symmetric and positive definite')
end

% initial conditions
%p=q;
p=zeros(mq,nq);
k=zeros(nb,n);
nrm=0;e=1;

while e>del
%  I=inv(r+b'*p*b)*b'*p;
   I=(r+b'*p*b)\(b'*p);
   k=I*a;
   p=a'*(p-p*b*I)*a+q;
   nrmp=nrm;
   nrm=norm(k);
   if nrm~=0, e=abs(nrmp-nrm)/nrm; else, e=1; end
end;

% end of m-file