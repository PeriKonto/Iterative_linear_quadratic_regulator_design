function T=Tmatrix(A_nm, A_min,at,bt) 

%This function calculates a T such as a transformation from Non-minimal Observable
%to minimal state vector would be possible 
%Author: Peri Kontoroupis
%date: 27/07/01
%Problems arose with the B values...
%B works only for two zeros in front...
%Must understand that we don't want the b1*z^-1
[I,J]=find(bt==0);

I=length(I);
if I >=2  
%disp('it goes through here');
bt=unpad(bt,0,'b');
%Observable form transformation
b=bt((I-1):length(bt));%truncated forms of the selected model, should check
else
bt=unpad(bt,0,'b');
b=bt((I+1):length(bt));
end
a=at(3:length(at));%in each case a^-1 should be discarted (1+a1*z^-1+a2*z^-2...)


[N1,M1]=size(A_nm);
[N2,M2]=size(A_min);
%M1=M1+1;
T=zeros(N2,M1);

T(1)=1;

ia=length(a);

ib=length(b);


for i=2:N2   %First line reserved with 0s go line by line until end

   l=i-1;
   %for l=1:ia    %find a GENERAL FORM including  both ia and ib
   %A values
   inter1=[-a(l:length(a))];
   ix=length(inter1);
   
          if ix<ia
             Val=ia-ix; %These times zeros are required
             %create a new array containing those zeros
             Zer=zeros(1,Val);
             %added it up to the original matirx
             inter1final=[inter1 Zer];
          else
             inter1final=inter1;
          end  
        
          inter1final;
          
  L=i-1; 
   
 %  for L=1:ib
      
   %L; try same for variable l
   inter2=[b(L:length(b))];
   iy=length(inter2);
    
         if ix<ib
            Val=ib-iy; %These times zeros are required
            %create a new array containing those zeros
            Zer=zeros(1,Val);
            %added it up to the original matirx
            inter2final=[inter2 Zer];
        else
            inter2final=inter2;
         end  
       
         inter2final;
         
Total=[inter1final inter2final];
T((l+1),:)=[0 Total];
   
  end   


