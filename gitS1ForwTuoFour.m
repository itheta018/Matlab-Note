clear;clc;close all;
clear all
close all
clc

%%%% 4 阶extended模型

A=[0 1 0 0;
    -28 -3.67 0 0;
    0 0 -1 pi/5;
    0 0 pi/5 -1];
Atest =  [-0.1 pi/5;
    pi/5 -0.1];
B1=[0;3.33; 0; 0];B2=[0;3.33; 0; 0];
% 
% 
% M=1396;
% Iz=1750;
% a=1.25;
% b=1.32;
% Cf=100700;
% Cr=86240;
% G=20.04;
% % vmax=20;
% % vmin=10;
% vx = 10;
% 
% a11=-(Cf+Cr)/(M*vx);
% a12=-(a*Cf-b*Cr)/(M*vx);
% a21=-(a*Cf-b*Cr)/(Iz*vx);
% a22=-(a^2*Cf+b^2*Cr)/(Iz*vx);
% a23 = (Cf+Cr)/(M);
% a43 = (a*Cf-b*Cr)/(Iz);
% 
% A=[0 1 0 0; 0 a11  a23  a12; 0 0 0 1; 0 a21 a43 a22];
% B1=[0;Cf/(M*G);0;a*Cf/(Iz*G)]; 
% B2=[0;Cf/(M*G);0;a*Cf/(Iz*G)]; 
Q1=diag([1,2,1,1]);R1=1;
Q2=diag([2,1,1,1]);R2=1;

P1=[1 0.1 1 1;
    0.1 1 1 1;
    1 1 1 1;
    1 1 1 1];
P2=[1 0.1 1 1;
    0.1 1 1 1;
    1 1 1 1;
    1 1 1 1];
ite=1;%iteMax=10;%20000;
iteMax=100;
maxN=1;
%while ((maxN > 10e-3 )&& (ite<iteMax)  ) ???
%while  maxN > 10e-2
while  ite<iteMax
A22=A-B1*inv(R1)*B1'*P1;
A11=A-B2*inv(R2)*B2'*P2;

[K2,P2,E2] = lqr(A22, B2, Q2, R2);
[K1,P1,E1] = lqr(A11, B1, Q1, R1);

A1L=A-B2*inv(R2)*B2'*P2;
A2L=A-B1*inv(R1)*B1'*P1;

N1L=A1L'*P1+P1*A1L-P1*B1*inv(R1)*B1'*P1+Q1;
N2L=A2L'*P2+P2*A2L-P2*B2*inv(R2)*B2'*P2+Q2;
maxN(ite) = max( norm(N1L,inf) ,norm(N2L,inf));
ite=ite+1
end
P1
P2
K1=inv(R1)*B1'*P1
K2=inv(R2)*B2'*P2
plot(maxN);

% K1=[0.393319893190329,0.273132129507922,6.846359012049240,0.370579686724889];
% K2=[1.074569931823545,0.173956392036126,6.040823516941583,0.355848388389530];

% A=[0 1 0; 0 0 1;-35 -27 -9];
% B=[0 0 1]';
% Q=[1 0 0;0 1 0;0 0 1];
% R=[1];
% 
% [K P E] = lqr(A, B, Q, R)
% 
%  eigenvalue=diag(P);
% 
%  lamda=max(eigenvalue)%P的最大特征值
% A=[0 1;0 0];
% 
% B1=[0;1];B2=[1;1];
% 
% R1=0.5;R2=2;
% Q1=[6 10;10 18];
% Q2=[6 13;13 30]

% %%%%%%%%%%%%%%%%% 2 order
% A=[0 1;-28 -3.67];B1=[0;3.33];B2=[0;3.33];
% R1=1;R2=1;
% Q1=[1 0;0 1];
% Q2=[2 0;0 1];
% 
% P1=[1 0.1;0.1 1];P2=[1 1;1 1];
% ite=1;%iteMax=10;%20000;
% iteMax=100;
% maxN=1;
% %while ((maxN > 10e-3 )&& (ite<iteMax)  ) ???
% %while  maxN > 10e-2
% while  ite<iteMax
% A22=A-B1*inv(R1)*B1'*P1;
% A11=A-B2*inv(R2)*B2'*P2;
% 
% [K2,P2,E2] = lqr(A22, B2, Q2, R2);
% [K1,P1,E1] = lqr(A11, B1, Q1, R1);
% 
% A1L=A-B2*inv(R2)*B2'*P2;
% A2L=A-B1*inv(R1)*B1'*P1;
% 
% N1L=A1L'*P1+P1*A1L-P1*B1*inv(R1)*B1'*P1+Q1;
% N2L=A2L'*P2+P2*A2L-P2*B2*inv(R2)*B2'*P2+Q2;
% maxN(ite) = max( norm(N1L,inf) ,norm(N2L,inf));
% ite=ite+1
% end
% P1
% P2
% K1=inv(R1)*B1'*P1
% K2=inv(R2)*B2'*P2
% plot(maxN);

%%%  K1 =    0.3429    0.2872  
%%%  K2 =    0.1698    0.8539

% K1 = [0.058445855837018,0.323523976928711];
% K2 = [0.117295148012668,0.335866825739659];
