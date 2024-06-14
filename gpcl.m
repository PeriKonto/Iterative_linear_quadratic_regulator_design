function [acl, bcl]=gpcl(a, b, S, R, W)
% GPCL  Closed loop transfer functions for GPC.
%
% [acl,bcl]=gpcl(a,b,s,r,w)
%
% b/a: System polynomials (trunciated form).
% r,s,w: GPC polynomials (see 'help gpc').
%
% bcl/acl: Command to output (Yd => y).

% James Taylor
% 09/08/1999

if nargin<5
  W=[];
end

bt=[0 b];
at=[1 a];

bcl=bt;
if ~isempty(W)  % command filter
  bcl=conv(bcl, W);
end
acl=polyadd(conv([1 -1], conv(R, at)), conv(S, bt));

% ensure denominator starts with unity
bcl=bcl/acl(1);
acl=acl/acl(1);

% end of m-file