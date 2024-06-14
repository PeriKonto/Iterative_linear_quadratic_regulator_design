function [F,g,d,h]=nmssform(a,b,inc);
% function [F,g,d,h]=nmssform(a,b);
% returns parameters of the NMSS description
% of a discrete system given by transfer function
% z polynomials a and b where a=1 + a1*z^-1 +..+an*z^-n
%
%  inc=1  forces the routine
%  to return the state space form without the integrated
%  control error state Z.

% Copyright (c) 1999 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor

% Revision history :
%   06/07/1999, JT, combined toolbox

ip=length(a);iq=length(b);
m=ip+iq;
% F matrix
z=ones(1,m-1);
F=diag(z,-1);
F(ip+1,:)=zeros(1,m);
z=[-a b(2:iq) 0];
F(1,:)=z;
F(m,:)=-z;
F(m,m)=1;
%
h=zeros(1,m);h(1)=1;
%
d=h';d(1)=0;d(m)=1;
%
g=h';g(1)=b(1);g(ip+1)=1;g(m)=-b(1);
%
if nargin==2, return;end
%if inc
  F=F(1:m-1,1:m-1);
  g=g(1:m-1);h=h(1:m-1);d=d(1:m-1);
%end
