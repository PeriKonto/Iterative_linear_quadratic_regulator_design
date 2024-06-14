%Program used for calculating root loci from a given transfer function.
%Three state space forms are used: namely observable,controllable forms, 
%as seen within the minimal state and the Lancaster development of non-minimal
%In conclusion all states are assumed available for measurement so there is no
%need using an observer at this point. Nevertheless future developments will
%take it uner consideration. Another future development will be the inclusion
%of an integral action within the NMSS form.
%Program includes the following functions:

%<-> tmod2 ........ Contains various choices of models to be used
%<-> dtf2ss .......  
%<-> nmssform .....
%<-> Tmatrix ...... 

%Based on previous work of J.C.Taylor 
%Program written by: Peri Kontoroupis
%----26-07-01---- System & Control Group
clc
clear all
close all
disp('-------Selected Model :-------')

[at, bt]=tmod2(7); %<----------------Choose your model from here!

b=bt(2:length(bt));%truncated forms of the selected model, should check
a=at(2:length(at));%whether really required...

%Some useful information concerning the model
sys=tf(bt,at,'variable','z^-1')%A last check for the selected model
disp('---poles---'); roots(at) %poles function gives an extra pole when system isn't coprime
disp('---zeros---'); zero(sys)

% Model Representation in Various State Space forms:
% (1)controllable canonical minimal state space form based on bt/at
[A_min, B_min, C_min]=dtf2ss(bt, at);
A_mincon=A_min; %Keep the states as is for use on plotting of Controllable form
B_mincon=B_min;
C_mincon=C_min;
D_mincon=zeros(size(C_mincon))';

% (2)observable canonical minimal state space form based on bt/at
A_min=A_min';   % transpose to convert to observable canonical form
B_min=zeros(size(A_min, 1), 1);
B_min(end-length(b)+1:end)=b';  % create observable canonical G matrix
C_min=zeros(size(B_min))'; % observation matrix
D_min=zeros(size(C_min))';
C_min(1)=1;

% (3)NMSS based on bt/at
[A_nm, B_nm, D_nm, C_nm]=nmssform(a, b, 1);

% Checking that all three state space representations are correct
disp('-------Controllable form----------');
A_mincon
B_mincon
C_mincon
disp('');
disp('--Contr. Poles Check--');
eig(A_mincon)
disp('');
disp('--------Observable form-----------');
A_min
B_min
C_min
disp('');
disp('--Observ. Poles Check--');
eig(A_min)
disp('');
disp('---------Non-minimal form---------');
A_nm
B_nm
C_nm
disp('');
disp('--Nonmin. Poles Check--');
eig(A_nm)
disp('');

%
% stochastic root loci for variations in Q_min, Q_mincon and Q_nm
% (Q_min:Observable, Q_mincon:Controllable, Q_nm:Nonminimal)

% cost function weightings for the observable
Q_min=C_min'*C_min;
R_min=0.001;

% cost function weightings
Q_nm=zeros(size(A_nm));
Q_nm(1:length(Q_min), 1:length(Q_min))=Q_min; %Provide Q_min with correct dimensions
R_nm=R_min;  % for equivalence with minimal

disp('------- Calculating stochastic root loci -------')

MC=100  % number of realisations <---Can be changed accordinally (1000 usual) 
Q_range=100  % weighting elements lie within zero to Q_range

% Loop done for Non-minimal since bigger, thus minimal forms would be selected
% from it. Future development the inclusion of a T transform. matrix

for ff=1:MC

  % progress
  if (MC>100) & (~rem(ff, 100)), disp(ff), end

  % random state weighting matrix
  if ff>1
    if 0  % diagonals only
       qd=rand(size(Q_nm, 1), 1)*Q_range+eps;  
       Q_nm=diag(qd); % eps: Floating point relative accuracy
    else
      QQ=rand(size(Q_nm))*Q_range+eps;
      QQ=triu(QQ);
      Q_nm=QQ'*QQ;  % ensures positive definite
      Qs{ff}=Q_nm;
    end
  end
  
%[i j]=size(Q_nm);  
%Diag_nm=diag(Q_nm);
%fact=0.01;
%Q_Q=fact*diag(Diag_nm);
%Q_Q(i,j)=Q_Q(i,j)*0.01;
%%Q_mincon=Q_min+Q_Q;  
%Q_nm=Q_Q;

% non-minimal (similar to above for each realisation)
K_nm=dlqri(A_nm, B_nm, Q_nm, R_nm, 1e-8);
F_closed = (A_nm-B_nm*K_nm);
d_closed=zeros(size(C_nm,1),size(D_nm,2));
sys = ss(F_closed,B_nm,C_nm,d_closed,-1);
r_nm=pole(sys);
c_nms(ff, 1)=length(find(abs(r_nm)<0.01));
r_nms(ff, :)=r_nm';

% minimal (Observable form)
Q_min=Q_nm(1:size(Q_min, 1), 1:size(Q_min, 1));
K_min=dlqri(A_min, B_min, Q_min, R_min, 1e-8);
F_closed = (A_min-B_min*K_min);
d_closed=zeros(size(C_min,1),size(D_min,2));
sys2 = ss(F_closed,B_min,C_min,d_closed,-1);
r_min=pole(sys2);
c_mins(ff, 1)=length(find(abs(r_min)<0.01));
r_mins(ff, :)=r_min';

T=Tmatrix(A_nm,A_min,at,bt);
%T matrix formulation: direct equivalence between NMSS and Minimal


% minimal (Controllable form)
Q_mincon=Q_min; %for equivalence to observable
K_mincon=dlqri(A_mincon, B_mincon, Q_mincon, R_min, 1e-8);
F_closed = (A_mincon-B_mincon*K_mincon);
d_closed=zeros(size(C_mincon,1),size(D_mincon,2));
sys3 = ss(F_closed,B_mincon,C_mincon,d_closed,-1);
r_mincon=pole(sys3);
c_mincon(ff, 1)=length(find(abs(r_mincon)<0.01));
r_minscon(ff, :)=r_mincon';

end


% plotting utility
disp('The T matrix (observable usage) is: ');
T
disp('');
disp('Provide a check using T matrix: ');
K_nm
K_check=K_min*T   %It wouldn't give the same result since K is varying (~Q) 
K_min      
fig('minimal observable form')
zgrid; axis('square')
plot(real(r_mins), imag(r_mins), 'k.')
fig('non-minimal form')
zgrid; axis('square')
plot(real(r_nms), imag(r_nms), 'k.')
fig('minimal controllable form')
zgrid; axis('square')
plot(real(r_minscon), imag(r_minscon), 'k.')
fig('all three together')
zgrid; axis('square')
plot(real(r_mins), imag(r_mins), 'r.')
plot(real(r_nms), imag(r_nms), 'k.')
plot(real(r_minscon), imag(r_minscon), 'b.')

