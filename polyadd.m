function c=polyadd(a,b,l)
% POLYADD  Adds two polynomials.
%
% c=polyadd(a,b)
%
% Adds polynomials a and b assuming that higher order
% coefficients are placed in last positions of a and b 
% (as in conv and deconv).

% James Taylor, 06/08/1999
% Wlodek Tych, 1990

ma=size(a);mb=size(b);
if (min(ma)>1)|(min(mb)>1)
  disp('Polyadd: data not polynomials !');
else

 if isempty(a), c=b;return;end
 if isempty(b), c=a;return;end

 if nargin==2  % poly. in z^-1
   mab=ma.*mb;       % check for row-column 
   ma=max(ma);mb=max(mb);
   if all(mab~=1), a=a'; end
   if ma==mb, c=a+b;return, end
   if ma>mb, c=a+[b zeros(1,ma-mb)]; return,end
   if ma<mb, c=b+[a zeros(1,mb-ma)]; return,end
 elseif nargin==3  % poly in z, s or delta
   mab=ma.*mb;       % check for row-column 
   ma=max(ma);mb=max(mb);
   if all(mab~=1), a=a'; end
   if ma==mb, c=a+b;return, end
   if ma>mb, c=a+[zeros(1,ma-mb) b]; return,end
   if ma<mb, c=b+[zeros(1,mb-ma) a]; return,end
 end
end

% end of m-file