             
             
clc
clear all

n=input('number of  unground nodes:')

% el=number of elements
% r=number of resistor
% c=number of capasitor
% L=number of inductor
% p=number of current source
% m=number of voltage source
% s=number of ccc
% b=number of ccv 
% a=number of vcc
% e=number of vcv
% ml=number of coupling inductor
            
             
E=load('input.txt')

el=length(find(E(:,1)))
r=length(find(E(:,1)==1))
c=length(find(E(:,1)==2))
L=length(find(E(:,1)==3))
p=length(find(E(:,1)==4))
m=length(find(E(:,1)==5))
f=length(find(E(:,1)==6))
b=length(find(E(:,1)==7))
a=length(find(E(:,1)==8))
e=length(find(E(:,1)==9))
ml=length(find(E(:,1)==10))

syms t
syms x
syms s
% modified node matrix
G=zeros(n,n);
P=0
for i=1:el
 if E(i,1)==1  
     
    
    for l=1:n
        P=0;
      for k=1:r
        if E(k,3)==l || E(k,4)==l
            if E(k,2)~=0
            P=1/E(k,2)+P;
            else
              1/E(k,2)==0;
            end
        end
      end 
       G(l,l)=P;
    end
    for k=1:r
          if E(k,3)~=0 && E(k,4)~=0
            if E(k,2)~=0
              G(E(k,3),E(k,4))=-1/E(k,2);
              G(E(k,4),E(k,3))=-1/E(k,2);
            else
              G(E(k,3),E(k,4))=0;
              G(E(k,4),E(k,3))=0; 
            end
          end
      end
 end
end
G
% capasitor
G1=zeros(n,n);
M1=0
syms  s
for i=1:el
 if E(i,1)==2  
     
    for l=1:n
        M1=0;
      for k=r+1:1:r+c
        if E(k,3)==l || E(k,4)==l
            
            M1=E(k,2)+M1;
           
            
        end
      end 
       G1(l,l)=M1;
    end
    for k=r+1:1:r+c
          if E(k,3)~=0 && E(k,4)~=0
              G1(E(k,3),E(k,4))=-E(k,2);
              G1(E(k,4),E(k,3))=-E(k,2);
            
          end
      end
 end
end
G1
G=G+(s*G1)


%VCC
G2=zeros(n,n);
for i=1:el
 if E(i,1)==8
   
    for j=1:a
        
        if E(i,3)~=0 &&  E(i,5)~=0  
           G2(E(i,5),E(i,3))=E(i,2);
        end
         
       if E(i,4)~=0 && E(i,6)~=0
          G2(E(i,6),E(i,4))=E(i,2);
         end
        
         if E(i,3)~=0 &&  E(i,6)~=0
            G2(E(i,6),E(i,3))=-E(i,2); 
         end
        
         if E(i,4)~=0 &&  E(i,5)~=0 
             G2(E(i,5),E(i,4))=-E(i,2);
         end
    end
 end
end
G=G+G2

% voltage source
B=zeros(n,m);
for i=1:el
 if E(i,1)==5
    
         if E(i,4)~=0
          B(E(i,4),i-(r+c+L+p))=-1;
         end
        
         if E(i,3)~=0
           B(E(i,3),i-(r+c+L+p))=1;
        end
    
 end
end
B

C = transpose(B)
D= zeros(m,m)
A=[G B;C D]


%Z contain the Is and Vs
%Z=[I;V]

% INDUCTOR
Q=zeros(n+m,L);
for i=1:el
 if E(i,1)==3
     
    if E(i,3)~=0
           Q(E(i,3),i-(r+c))=1;
    end
        
    if E(i,4)~=0
          Q(E(i,4),i-(r+c))=-1;
    end
        
         
        
 end
end
Q
R = transpose(Q)
S= zeros(L,L)
for i=1:el
 if E(i,1)==3
   for j=1:L
  S(j,j)=-E(i,2)

   end
 end
end
S=s*S

% A=[A Q;R S]


% COUPLING INDUCTOR
W=zeros(n+m+L,(2*ml));
for i=1:el
 if E(i,1)==10
     
    if E(i,3)~=0
           W(E(i,3),i-(r+c+L+p+m+f+b+a+e))=1;
    end
        
    if E(i,4)~=0
          W(E(i,4),i-(r+c+L+p+m+f+b+a+e))=-1;
    end
    
    if E(i,5)~=0
           W(E(i,5),i-(r+c+L+p+m+f+b+a+e))=1;
    end
        
    if E(i,6)~=0
          W(E(i,6),i-(r+c+L+p+m+f+b+a+e))=-1;
    end
          
 end
end
W
X = transpose(W)
Y= zeros((2*ml),(2*ml))
for i=1:el
 if E(i,1)==10
   for j=1:ml
  Y(j,j)=-E(i,2)
  Y(j+1,j+1)=-E(i,8)
  Y(j+1,j)=-E(i,9)
  Y(j,j+1)=-E(i,9)
   end
 end
end
Y=s*Y
% A=[A W;X Y]

% CURRENT SOURCE
I(n,1)=vpa(0)
for i=1:el
 if E(i,1)==4
     if E(i,11)==0
         if E(i,3)~=0
          I(E(i,3),1)=-E(i,2);
         end
          if E(i,4)~=0
          I(E(i,4),1)=E(i,2);
          end
     else if E(i,11)==1
           if E(i,3)~=0
          I(E(i,3),1)=-laplace(E(i,2)*x^0);
           end
          if E(i,4)~=0
          I(E(i,4),1)=laplace(E(i,2)*x^0);
          end
         else if E(i,11)==2
         if E(i,3)~=0
          I(E(i,3),1)=-laplace(sin(E(i,2)*t));
         end
          if E(i,4)~=0
          I(E(i,4),1)=laplace(sin(E(i,2)*t));
          end
             end
         end
     end
    end
end
I
% initial states for capasitor
I1=zeros(n,1)
for i=1:el
 if E(i,1)==2
   
         if E(i,3)~=0
          I1(E(i,3),1)=E(i,2)*E(i,7);
         end
          if E(i,4)~=0
          I1(E(i,4),1)=-E(i,2)*E(i,7);
          end
          
    end
end
I1
I2=I+I1

%Z=zeros(L,1)

%INITIAL STATE FOR INDUCTOR
I0=sym(zeros(L,1))
for i=1:el
 if E(i,1)==3
 
          I0(i-(r+c),1)=-E(i,2)*E(i,7);
end
end 
% INITIAL STATE FOR COUPLING INDUCTOR

% VOLTAGE SOURCE
V=sym(zeros(m,1))
for i=1:el
 if E(i,1)==5
      if E(i,11)==0
          V(i-(r+c+L+p),1)=E(i,2)
      elseif E(i,11)==1
           V(i-(r+c+L+p),1)=laplace(E(i,2)*x^0)
          elseif E(i,11)==2
         V(i-(r+c+L+p),1)=laplace(sin(E(i,2)*t))
      end
 end
end

V   

 
% Z=[I2;I0;V]
for i=1:el
if E(i,1)==2
    Z=I2
end
end

for i=1:el
if E(i,1)==3
    A=[A Q;R S]
    Z=[I2;I0]
end
end 

for i=1:el
if E(i,1)==5
    Z=[I2;I0;V]
end
end

for i=1:el
 if E(i,1)==10
    A=[A W;X Y]
    I3=zeros((2*ml),1)
    Z3=[I2;I0;V;I3]
    Z=Z3
 end
end

% X=inv(A)*Z

%VCV
M=zeros(n+m+L+(2*ml),e);
for i=1:el
 if E(i,1)==9
       for j=1:e
         if E(i,4)~=0
          M(E(i,4),j)=-1;
         end
        
         if E(i,3)~=0
           M(E(i,3),j)=1;
         end
       end
    end
end
N=transpose(M)
N1=zeros(e,n+m+L+(2*ml));
for i=1:el
 if E(i,1)==9
   for j=1:e
         if E(i,5)~=0
          N1(j,E(i,5))=E(i,2);
         end
        
         if E(i,6)~=0
           N1(j,E(i,6))=-E(i,2);
         end
   end
    end
end
N=N+N1
O=zeros(e,e)

%A=[G B;C D]
for i=1:el
if E(i,1)==9
    A=[A M;N O]
    L3=zeros(e,1)
    Z=[Z;L3]
end
end
%X=inv(A)*Z

%CCC
MM=zeros(n+m+L+(2*ml)+e,f);
for i=1:el
 if E(i,1)==6
    
       for j=1:f
         if E(i,4)~=0
          MM(E(i,4),j)=-E(i,2);
         end
        
         if E(i,3)~=0
           MM(E(i,3),j)=E(i,2);
         end
       end
    end
end

MM1=zeros(n+m+L+(2*ml)+e,f);
for i=1:el
 if E(i,1)==6
   for j=1:f
         if E(i,5)~=0
          MM1(E(i,5),j)=1;
         end
        
         if E(i,6)~=0
           MM1(E(i,6),j)=-1;
         end
   end
 end
end
MM=MM+MM1
NN=transpose(MM1)
OO=zeros(f,f)
for i=1:el
 if E(i,1)==6
    A=[A MM;NN OO]
    L1=zeros(f,1)
    Z=[Z;L1]
 end
end

%X=inv(A)*Z

%CCV
MMM=zeros(n+m+L+(2*ml)+e+f,b+1);
for i=1:el
 if E(i,1)==7
    
       for j=1:b
           if E(i,3)~=0
            MMM(E(i,3),j+1)=1;
           end
           if E(i,4)~=0
            MMM(E(i,4),j+1)=-1;
           end
           if E(i,5)~=0
            MMM(E(i,5),j)=1;
           end
           if E(i,6)~=0
            MMM(E(i,6),j)=-1;
           end
        end
         
   end
end  
NNN=transpose(MMM)
OOO=zeros(b+1,b+1)
for i=1:el
 if E(i,1)==7
   for j=1:b
  OOO(j+1,j)=-E(i,2)

   end
 end
end
for i=1:el 
 if E(i,1)==7
    A=[A MMM;NNN OOO]
    L2=zeros(b+1,1)
    Z=[Z;L2]
 end
end
 
X=inv(A)*Z
Z
        

% SPECIFY  H(s) AND SKETCH H(jw) AND h(t)           
  for i=1:el
      if E(i,4)~=0
                 if E(i,10)==1
                     
                      if E(i,1)==1
                         I4=(X(E(i,3),1)-X(E(i,4),1))/(E(i,2))
                      elseif E(i,1)==2
                         I4=(X(E(i,3),1)-X(E(i,4),1))*s*E(i,2);
                      else  E(i,1)==3
                         I4=(X(E(i,3),1)-X(E(i,4),1))/(s*E(i,2));
                      end
                         H1=I4;
                     
                     elseif E(i,10)==2
                         H2= E(i,2)
                         
                         
                     elseif E(i,10)==3
                     V1=X(E(i,3),1)-X(E(i,4),1);
                     H1=V1;
                     
                     elseif E(i,10)==4
                     H2=E(i,2);
                     
                     elseif E(i,10)==5
                          if E(i,1)==1
                         I4=(X(E(i,3),1)-X(E(i,4),1))/E(i,2);
                     elseif E(i,1)==2
                         I4=(X(E(i,3),1)-X(E(i,4),1))*(s*E(i,2));
                     else  E(i,1)==3
                         I4=(X(E(i,3),1)-X(E(i,4),1))/(s*E(i,2));
                          end
                     H1=I4;
                      
                     elseif E(i,10)==6
                          H2=E(i,2);
                          
                          
                     elseif E(i,10)==7
                     V1=X(E(i,3),1)-X(E(i,4),1);
                     H1=V1;
                    
                     elseif E(i,10)==8
                         H2=E(i,2);
                     
                 end
       elseif E(i,4)==0
                    if E(i,10)==1 
                     
                        if E(i,1)==1
                         I4=X(E(i,3),1)/E(i,2)
                      elseif E(i,1)==2
                         I4=(X(E(i,3),1))*s*E(i,2);
                      else  E(i,1)==3
                         I4=(X(E(i,3),1))/(s*E(i,2));
                        end
                         H1=I4
                     
                     elseif E(i,10)==2
                         H2= E(i,2)
                         
                         
                 elseif E(i,10)==3
                     V1=X(E(i,3),1);
                     H1=V1;
                 
                     elseif E(i,10)==4
                     H2=E(i,2);
                     
                     elseif E(i,10)==5
                          if E(i,1)==1
                         I4=(X(E(i,3),1))/E(i,2);
                     elseif E(i,1)==2
                         I4=(X(E(i,3),1))*(s*E(i,2));
                     else  E(i,1)==3
                         I4=(X(E(i,3),1))/(s*E(i,2));
                          end
                     H1=I4;
                     
                     elseif E(i,10)==6
                          H2=E(i,2);
                          
                     elseif E(i,10)==7
                           V1=X(E(i,3),1);
                     H1=V1;
                    
                     elseif E(i,10)==8
                         H2=E(i,2);
                    end
                 end
      end

  
         
           
   
 syms answer1  answer2 yes t            
             
rsp=input('which response do you want?')

if rsp==1
%ZERO STATE RESPONSE
I0=zeros(L,1);
I1=zeros(n,1);
I2=I+I1
I3=zeros((2*ml),1)
Z3=[I2;I0;V;I3]
Xzeros=inv(A)*Z
end


if rsp==2
%ZERO INPUT RESPONSE
I=zeros(n,1);
V=zeros(m,1);
I2=I+I1
I3=zeros((2*ml),1)
Z3=[I2;I0;V;I3]
Xzeroi=inv(A)*Z
end



% FIND NATURAL FREQUENCIES 
syms s j w
   if rsp==2
       det(G)
      NF=solve(det(G))
      NFB=vpa(NF)
   end   
   
% SPECIFY  H(s) AND SKETCH H(jw) AND h(t) continiue
              if rsp==1
                answer1=input('Do you want H(s)?')
             end 
  syms s j w           
             if answer1==yes
             H=H1/H2
             ezplot(H);
             hold on;
             h=ilaplace(H,t)
             ezplot(h);
             title 'impulse response'
             xlabel('t axis')
             ylabel('h axis')
             
             end
             
% determine the current and voltage for each element 
answer2=input('DO you want voltage and current response?')
             i3=1;
             if answer2==yes
             for i3=1:el
                 if E(i,4)~=0
                   if E(i3,12)==1
                      if E(i3,1)==1
                         II=(X(E(i,3),1)-X(E(i,4),1))/(E(i,2));
                     elseif E(i3,1)==2
                         II=(X(E(i,3),1)-X(E(i,4),1))*s*E(i,2);
                     else  E(i3,1)==3
                         II=(X(E(i,3),1)-X(E(i,4),1))/(s*E(i,2));
                      end
                      ii=ilaplace(II)
                   elseif E(i3,12)==2
                          VV=X(E(i3,3),1)-X(E(i3,4),1);
                          vv=ilaplace(VV)
                
                   end
                   elseif E(i,4)==0
                       if E(i3,12)==1
                      if E(i3,1)==1
                         II=(X(E(i,3),1))/E(i,2);
                     elseif E(i3,1)==2
                         II=(X(E(i,3),1))*s*E(i,2);
                     else  E(i3,1)==3
                         II=(X(E(i,3),1))/(s*E(i,2));
                      end
                      ii=ilaplace(II)
                       elseif E(i3,12)==2
                          VV=X(E(i,3),1);
                          vv=ilaplace(VV)
                
                   end
                 end
             end
             end
           
            
