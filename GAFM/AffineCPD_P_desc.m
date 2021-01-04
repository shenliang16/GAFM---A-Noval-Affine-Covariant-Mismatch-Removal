function [P1,P,E]=AffineCPD_P_desc(Xi,Ti, sigma ,omiga, D_desc, sigma2, Sigpower)

[N, D]=size(Xi); 
X1=Xi(:,1:2);   X2=Xi(:,3:4);   X3=Xi(:,5:6);   %X4=des1;
T1=Ti(:,1:2);   T2=Ti(:,3:4);   T3=Ti(:,5:6);   %T4=des2;

son2= ((2*pi)^(D/2 + Sigpower)) * (sigma * sigma * sigma * sigma2.^Sigpower) * (omiga/(1-omiga)) / 3^3;

sontemp = (X1-T1).*(X1-T1)/(2*sigma)+...
          (X2-T2).*(X2-T2)/(2*sigma)+...
          (X3-T3).*(X3-T3)/(2*sigma); 
D_desc_sigma = D_desc/(2*sigma2);

sontemp1 = sum(sontemp,2); 
son1 = exp(-sontemp1-D_desc_sigma);  % ∏≈¬ œÚ¡ø
son = son1+son2;
P1 = son1./son;
E = sum(P1.*sontemp1);

P=diag(P1); Np=sum(P1);  
E = E + 3*Np*log(sigma) + Sigpower*Np*log(sigma2);
E = E - Np*log(1-omiga) - (N-Np)*log(omiga);   

end