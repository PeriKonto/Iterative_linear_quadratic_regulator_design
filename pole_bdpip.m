%pole_bdpip
%Funtion for computing the SVF gains k 
%--------------------------(2)---------------------------------
%-----------------Pole assignment design-----------------------
%-----------------Block Diagram analysis-----------------------
%k=(MS^T)^-1(d-p)
%d & p are vectors of coeffients for the desired closed/open loop char.poly
%S, M control.mtrix and lower triangular respectively
%req. system must be controllable
%by Peri Kontoroupis 21/01/02


function k=pole_bdpip(num,den,d_poles)
error(nargchk(3,3,nargin));

clc
clear all
close all

%for cases m <= n
%EXAMPLE (1)						%WORKING
%d_poles=[0.6 0.5 0.3 0.2]; %Rogers thesis p.36
%num=[0 1 0.5];
%den=[1 -1.2 0.35];

%EXAMPLE (2)						%WORKING
%d_poles=[.5 .5 .5 .5]; 		% n > m
%num=[0 0 -1 2];
%den=[1 -1.7 1];

%EXAMPLE (3)						%WORKING
%d_poles=[0.8 0.8]; 
%num=[0 2.25];
%den=[1 -0.713];
%--------------------------------------------------

sys=tf(num,den,'variable','z^-1') %system check

b=num(2:length(num)); a=den(2:length(den));
[F,g]=nmssform(a,b);

m=length(num); n=length(den); %p computation
