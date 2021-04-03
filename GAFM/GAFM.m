function [C, sigma, iter, Ti] = GAFM(X,Y, max_it, tol, viz, omiga, SiftRatio, P_thr)
% X, Y:         Size: N¡Á6;              Meaning: Prematched affine correspondence (using ASIFT, Hessian-Affine etc.),   
% max_it:       Maximum iteration
% tol:          tolerance,               1e-3 - 1e-5;
% viz:          visualization flag,      0/1
% omiga:        outlier ratio
% SiftRatio:    Size: N¡Á1;              Sift Ratio Values corresponding to {X, Y}
% P_thr:        Decision threshold.

X = double(X(:,[3,4,5,6,1,2]));   Y = double(Y(:,[3,4,5,6,1,2]));  

[X,Y]=LAF_to_pt3x3(X, Y);
[X,Y,normal]=AffineCPD_normalize(X,Y);
[C, ~, ~, sigma, iter, Ti]=our_affine2(X,Y, max_it, tol, viz, omiga, SiftRatio, P_thr);
end

function [C, A, t, sigma, iter, Ti]=our_affine2(X,Y, max_it, tol, viz, omiga, SiftRatio, P_thr)

[N, ~]=size(X); [M, ~]=size(Y); D=2;   
    
lamnda=0;     Sigpower=0.5;

D_desc=SiftRatio;       
D_desc = (D_desc-min(D_desc))/(max(D_desc)-min(D_desc));


%% Initialize sigma and Y
X1=X(:,1:2);   X2=X(:,3:4);   X3=X(:,5:6);
Y1=Y(:,1:2);   Y2=Y(:,3:4);   Y3=Y(:,5:6);
sigma=sum(sum((X-Y).^2))/(N*6); if sigma==0; sigma=1; end
sigmaDesc=0.02;%sigmaDesc=0.02;  %Sigpower=0.5;
Xi=X; Yi=Y; 
Ti=Yi;

% Optimization
iter=0; ntol=tol+10; L=1; %omiga=0.5;
while (iter<max_it) && (ntol > tol) && (sigma > 10*eps)
% for iter=1:200
    L_old=L;
%     [P1,Pt1, P,L]=cpd_P_ours(Xi,Ti, [sigma,sigma,sigma] ,omiga); st='';
    
    [P1,~,L]=AffineCPD_P_desc(Xi,Ti, sigma ,omiga, D_desc, sigmaDesc, Sigpower); st='';
    P1 = max(P1, 1e-8);
    Pt1=P1'; 
    P=diag(P1);
    
    ntol=abs((L-L_old)/L);
%     disp([' CPD Affine ' st ' : dL= ' num2str(ntol) ', iter= ' num2str(iter) ' sigma2= ' num2str(sigma)]);
  
    % Precompute 
    Np=sum(P1);
    
    mu_x=Pt1*(X1+X2+X3)/(3*Np);
    mu_y=(Y1+Y2+Y3)'*P1/(3*Np); mu_y=mu_y';
    
    Xm1=(X1-mu_x);     Xm2=(X2-mu_x);      Xm3=(X3-mu_x);
    Ym1=(Y1-mu_y);     Ym2=(Y2-mu_y);      Ym3=(Y3-mu_y);
    % Solve for parameters
    A1=Xm1'*P'*Ym1+Xm2'*P'*Ym2+Xm3'*P'*Ym3;  %PX
    A2=Ym1'.*repmat(P1',[2,1])*Ym1 + Ym2'.*repmat(P1',[2,1])*Ym2 + Ym3'.*repmat(P1',[2,1])*Ym3;  
    A=A1/A2; % B= B1 * inv(B2);
    t=mu_x'-A*mu_y';  
    
   %% ¸üÐÂsigma
    sigmaDesc = 2*sum(P1.*D_desc)/(2*Sigpower*Np);
    sigma1=(trace(Xm1'.*repmat(Pt1,[2,1])*Xm1)-2*trace(Xm1'*P'*Ym1*A')+trace(A*Ym1'.*repmat(P1',[2,1])*Ym1*A'))/Np;  
    sigma2=(trace(Xm2'.*repmat(Pt1,[2,1])*Xm2)-2*trace(Xm2'*P'*Ym2*A')+trace(A*Ym2'.*repmat(P1',[2,1])*Ym2*A'))/Np;  
    sigma3=(trace(Xm3'.*repmat(Pt1,[2,1])*Xm3)-2*trace(Xm3'*P'*Ym3*A')+trace(A*Ym3'.*repmat(P1',[2,1])*Ym3*A'))/Np;
    sigma=(sigma1+sigma2+sigma3)/6;

  
    Ti(:,1:2)=Y(:,1:2)*A'+repmat(t',[M 1]);
    Ti(:,3:4)=Y(:,3:4)*A'+repmat(t',[M 1]);
    Ti(:,5:6)=Y(:,5:6)*A'+repmat(t',[M 1]);
    
     iter=iter+1;
    if ~viz, cpd_plot_iter(X(:,5:6), Ti(:,5:6)); end;
    
end

      [P1,~,L]=AffineCPD_P_desc(Xi,Ti, sigma ,omiga, D_desc, 10, Sigpower); st='';
      indX=find(P1>P_thr);
      C=indX;
end
      
      