function [at, bt, a, b]=tmod(model)
% TMOD  Various discrete time transfer function models.
%
% [at,bt,a,b]=tmod(model)
%
% Model: model number 1-5, see below.
%
% bt/at: Polynomials in full form.
% b/a: Polynomials in trunciated form,
%        where bt=[0 b] and at=[1 a]
%
% model 1 : 1st order stable, 2 samples time delay.
%
% model 2 : 2nd order, non minimum phase, self oscillatory.
%
% model 3 : highly unstable 3rd order
%           based on continuous time model with the
%           following type of structure: (s-0.1)(s^2 - s + 1).
%             [bt,at]=ord2(2, -0.5);
%             at=conv(at, [1 -0.1]);
%             [bt,at]=c2dm(bt, at, 0.2);
%
% model 4 : non minimum phase 1st order stable
%             from appendix B, Clarke et al. (1987)
%             Automatica, 23, 137-148.
%
% model 5 : 1st order glasshouse model, unit time delay.
%
% model 6 : 1st order stable
%             from p32, R.R. Bitmead et al. (1990)
%             Adaptive Optimal Control, Prentice Hall.

% James Taylor
% 09/08/1999

if model<=1
  at=[1 -0.9];
  bt=[0 0 0.005];
elseif model<=2
  at=[1 -1.7 1];
  bt=[0 0 -1 2];
elseif model<=3
  at=[1 -3.3179 3.8359 -1.5220];
  bt=[0 0.0015 0.0065 0.0018];
elseif model<=4
  at=[1 -0.9];
  bt=[0 1 2];
elseif model<=5
  at=[1 -0.8738];
  bt=[0 2];
elseif model<=6
  at=[1 -0.7];
  bt=[0 0.9 -0.6];
end

a=at(2:length(at));
b=bt(2:length(bt));

% end of m-file