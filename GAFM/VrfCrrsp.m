function [Crr_GT,d1,d2] = VrfCrrsp(B_opt,  X1, X2, threshold)

    T=B_opt([2,1,3],[2,1,3]); 
    d2=[]; 
    threshold=max(5,threshold);
    [Crr_GT,d1]=DetGTusingT(X1,X2,threshold,T);
end
    
function [GrdTrth, Dis]=DetGTusingT(X,Y,thr,T)

    N=size(X,1);    X3=[X(:,1:2),ones(N,1)];     
    XT=(T*X3')';        TransformedX=XT(:,1:2)./XT(:,3);

    Dis=sqrt(sum((Y(:,1:2)-TransformedX).^2,2));

    GrdTrth=(Dis<thr);
    
end