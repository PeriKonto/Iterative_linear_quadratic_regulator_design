% Stochastic root loci for variations in the state weighting matrix
% comparing minimal and non-minimal state space designs
%
% This is the example used in:
%   Taylor, C.J., Chotai, A. and Young, P.C., The unification of state
%   space control system design by non-minimal state variable feedback, 
%   in preparation for submission to IJC (September 1999).

% James Taylor
% 02/09/1999
% Revised/Modifications by:
% Peri Kontoroupis
% 20/07/2001
clc
clear all
close all

%
% control model
%

% deterministic model
disp('------- Model -------')
[at, bt]=tmod2(8); %<----------------Choose your model from here!
%at=conv(at, [1 -0.9]);  %Model initially chosen is Different from original seen in Paper.
%bt=conv(bt, [1 -0.7]);
b=bt(2:length(bt)); %this is done to have correct forms in NMSS
a=at(2:length(at)); %SHOULD CHECK WHETHER USED ANYWHERE ELSE

sys=tf(bt,at,'variable','z^-1');
sys
disp('Or :')
tf(bt,at,'variable','z')
disp('---poles---');
%pole(sys)
roots(at)
disp('---zeros---');
zero(sys)
%roots(bt)
% stochastic model
%adt=conv(at, [1 -1]);
%ad=adt(2:length(adt));
%ad=at(2:length(at)); 
ad=a;
%bt %<-------------------
%at %<------------------- These do give the 4rth order model seen in paper by J.T
%adt

%
% regulator minimal LQ
%

% observable canonical minimal state space form based on bt/adt
%[A_min, B_min, C_min]=dtf2ss(bt, adt);  % returns controllable form
%-------------------------------------------
%%Added code by Peri Kontoroupis 20-08-01 %note 'at' changed to 'adt' as above
[A_min, B_min, C_min]=dtf2ss(bt, at);
%disp(' '); disp(' ');
%disp('----Model Chosen to carry the analysis---');
%sys2=tf(bt,at,'variable','z^-1');
%sys2%<--------------Added by PERI
%disp('---poles---');
%pole(sys2)
%roots(at)
%disp('--zeros--')
%disp(' ');
%zero(sys2)
%roots(bt)
%disp('');
%disp('Or :')
%tf(bt,at,'variable','z')
A_mincon=A_min; %Keep the states as is for use on plotting of Controllable form
B_mincon=B_min;
C_mincon=C_min;
D_mincon=zeros(size(C_mincon))';

disp('-------Controllable form----------')
A_mincon
B_mincon
C_mincon
%%end of addition---------------------------
disp('--Check--:');
eig(A_mincon)
A_min=A_min';  % transpose to convert to observable canonical form
B_min=zeros(size(A_min, 1), 1);
B_min(end-length(b)+1:end)=b';  % create observable canonical G matrix
C_min=zeros(size(B_min))';  % observation matrix
C_min(1)=1;
disp('--------Observable form-----------');
A_min
B_min
C_min
disp('--Check--:');
eig(A_min)

% cost function weightings
Q_min=C_min'*C_min;
R_min=0.01; %<--------------- Change this value
disp('');

% infinite recursion
disp('------- Infinite minimal LQ with deadbeat observer -------')
[K_min, P_min]=dlqri(A_min, B_min, Q_min, R_min, 1e-8);
%K_min
K_min_bitmead=-K_min;  % reverse sign for equivalence with Bitmead

% observer
M=C_min';  % Bitmead's deadbeat solution

% observer matrices
A_obs=(eye(size(M*C_min))-M*C_min)*A_min;
B_obs=(eye(size(M*C_min))-M*C_min)*B_min;

% convert bt/adt based control law from state space to
% transfer function, where num/den=C*inv(sI-A)*B+D
A_ss2tf=A_min-M*C_min*A_min+(B_min-M*C_min*B_min)*K_min_bitmead;
B_ss2tf=M;
C_ss2tf=K_min_bitmead;
D_ss2tf=zeros(size(C_ss2tf, 1), size(B_ss2tf, 2));
% use s-operator based ss2tf
[S_min, U_min]=ss2tf(A_ss2tf, B_ss2tf, C_ss2tf, D_ss2tf);

% convert to GPC structure
S_min=S_min(2:length(S_min));  % convert to backward shift
S_min=-S_min  % convert to negative feedback convention
U_min=unpad(U_min, 0, 'e')  % remove trailing zeros

% poles
[acl,bcl]=gpcl(a, b, S_min, U_min);
r_min=roots(acl) %is r_min used anywhere else ?
disp('---Bitmead Observer model---');
tf(U_min,S_min,'variable','z^-1')
disp('');
disp('---Actual model---');
tf(b,a,'variable','z^-1')
disp('')
%
% regulator non-minimal LQ based on differenced inputs
%

% regulator NMSS based on bt/adt
%[A_nm, B_nm, tmp, C_nm]=nmssform(ad, b, 1);
[A_nm, B_nm, tmp, C_nm]=nmssform(a, b, 1); 
disp('---Model used for NMSS (only for the correct State-Space)---')
sys3=tf(b,a,'variable','z^-1');
sys3
disp('---poles---');
%pole(sys3)
roots(a)
disp('------Non-minimal form--------');
A_nm
B_nm
C_nm
% cost function weightings
Q_nm=zeros(size(A_nm));
Q_nm(1:length(Q_min), 1:length(Q_min))=Q_min;
R_nm=R_min;  % for equivalence with minimal

% infinite recursion
disp('------- Infinite NMSS LQ (differenced inputs) -------')
P0_nm=Q_nm;  % initial P matrix
[K_nm, P_nm]=dlqri(A_nm, B_nm, Q_nm, R_nm, 1e-8);

% convert to GPC structure
S_nm=K_nm(1:length(ad))
U_nm=[1 K_nm(length(ad)+1:end)]

% poles
[acl,bcl]=gpcl(a, b, S_nm, U_nm); %<----------- check this
r_nm=roots(acl);%<----------- check this


%
% stochastic root loci for variations in Q_min and Q_nm
%

disp('------- Calculating stochastic root loci -------')

MC=1000  % number of realisations
Q_range=1000  % weighting elements lie within zero to Q_range

for ff=1:MC

  % progress
  if (MC>100) & (~rem(ff, 100)), disp(ff), end

  % random state weighting matrix
  if ff>1
    if 0  % diagonals only
      qd=rand(size(Q_nm, 1), 1)*Q_range+eps;
      Q_nm=diag(qd);
     
    else
      QQ=rand(size(Q_nm))*Q_range+eps;
      QQ=triu(QQ);
      Q_nm=QQ'*QQ;  % ensures positive definite
      Qs{ff}=Q_nm;
    end
  end

% non-minimal (similar to above for each realisation)
  %Added by PERI
  %Diag_nm=diag(Q_nm);
  %fact=0.1;
  %Q_Q=fact*diag(Diag_nm);
  %Q_Q(3,3)=Q_Q(3,3)*0.1;
  %%Q_mincon=Q_min+Q_Q;  
  %Q_nm=Q_Q;
  %--------------

  K_nm=dlqri(A_nm, B_nm, Q_nm, R_nm, 1e-8);
  S_nm=K_nm(1:length(ad));
  U_nm=[1 K_nm(length(ad)+1:end)];
  acl=gpcl(a, b, S_nm, U_nm);
  r=roots(acl);
  c_nms(ff, 1)=length(find(abs(r)<0.01));  
  r_nms(ff, :)=r';

  % minimal (similar to above for each realisation)
  Q_min=Q_nm(1:size(Q_min, 1), 1:size(Q_min, 1));
  %Diag_min=diag(Q_min);
  %fact=0.1;
  %Q_Q=fact*diag(Diag_min);
  %Q_Q(4,4)=Q_Q(4,4)*0.001;
  %Q_mincon=Q_min+Q_Q;  
  %Q_min=Q_Q;

  K_min=dlqri(A_min, B_min, Q_min, R_min, 1e-8);
  K_min_bitmead=-K_min;  % reverse sign for equivalence with Bitmead
  A_ss2tf=A_min-M*C_min*A_min+(B_min-M*C_min*B_min)*K_min_bitmead;
  C_ss2tf=K_min_bitmead;
  [S_min, U_min]=ss2tf(A_ss2tf, B_ss2tf, C_ss2tf, D_ss2tf);
  S_min=S_min(2:length(S_min));  % convert to backward shift
  S_min=-S_min;  % convert to negative feedback convention
  U_min=unpad(U_min, 0, 'e');  % remove trailing zeros
  acl=gpcl(a, b, S_min, U_min);
  r=roots(acl);
  c_mins(ff, 1)=length(find(abs(r)<0.01));  % allows for numerical ...
  r_mins(ff, :)=r';
  
  %-------------------------------------------
  %%Added code by Peri Kontoroupis 20-08-01
  % minimal controllable form
  %%Q_mincon=Q_nm(1:size(Q_min, 1), 1:size(Q_min, 1));
  %%Diag_min=diag(Q_min);
  %%fact=0.1;
  %%Q_Q=fact*diag(Diag_min);
  %%Q_Q(4,4)=Q_Q(4,4)*0.001;
  %%Q_mincon=Q_min+Q_Q;  
  %%Q_mincon=Q_Q;
  
  Q_mincon=Q_min;
  K_mincon=dlqri(A_mincon, B_mincon, Q_mincon, R_min, 1e-8);
  F_closed = (A_mincon-B_mincon*K_mincon);
  d_closed=zeros(size(C_mincon,1),size(D_mincon,2));
  sys = ss(F_closed,B_mincon,C_mincon,d_closed,-1);
  r_mincon=pole(sys);
  c_mins(ff, 1)=length(find(abs(r_mincon)<0.01));
  r_minscon(ff, :)=r_mincon';
  %%end of addition---------------------------
  
end

% plot separately
fig('minimal observable form GPC equivalent')
zgrid; axis('square')
plot(real(r_mins), imag(r_mins), 'k.')
fig('non-miminal form')
zgrid; axis('square')
plot(real(r_nms), imag(r_nms), 'k.')
%fig('both')
%zgrid; axis('square')
%plot(real(r_mins), imag(r_mins), 'r.')
%plot(real(r_nms), imag(r_nms), 'k.')

%-------------------------------------------
%%Added code by Peri Kontoroupis 20-08-01
fig('minimal controllable form')
zgrid; axis('square')
plot(real(r_minscon), imag(r_minscon), 'k.')
fig('all three together')
zgrid; axis('square')
plot(real(r_mins), imag(r_mins), 'r.')
plot(real(r_nms), imag(r_nms), 'k.')
plot(real(r_minscon), imag(r_minscon), 'b.')

%%end of addition---------------------------

%end of m-file