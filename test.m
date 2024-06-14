clc
clear all
close all
%(1)
   %at=[1 -0.7];
   %bt=[0 0.9 -0.6];
%(2)
   bt=[0 0 -1 2];
   at=[1 -1.7 1];
% (3)
  % bt=[0 2.25]; 
  % at=[1 -0.713];
% (4)  
  % at=[1 -2.6 2.53 -0.9];
  % bt=[0 0 -1 2.7 -1.4];
% (5)  

   % at=[1 -0.8 0.15];
   % bt=[0 0.5 -0.4]; 
 % (6) 
 
 	% at=[1 -1.2 0.35];;
   % bt=[0 1 0.5];
   
  sys=tf(bt,at,'variable','z^-1')  
  at=at(2:end);
  bt=bt(2:end);
 
  [F,g]=nmssform(at,bt)