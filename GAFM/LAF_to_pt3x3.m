function [Xp,Yp] = LAF_to_pt3x3(Xf, Yf)

if isempty(Xf)
   Xp = []; Yp=[]; 
   return; 
end

temp = sum(Xf);
if temp(3)>0.05*temp(5) 
    Xp = Xf;
    Yp = Yf;  
    return;
end
affptx(:).x=Xf(:,5)';           affpty(:).x=Yf(:,5)';
affptx(:).y=Xf(:,6)';           affpty(:).y=Yf(:,6)';
affptx(:).a11=Xf(:,1)';         affpty(:).a11=Yf(:,1)';
affptx(:).a12=Xf(:,3)';         affpty(:).a12=Yf(:,3)';
affptx(:).a21=Xf(:,2)';         affpty(:).a21=Yf(:,2)';
affptx(:).a22=Xf(:,4)';         affpty(:).a22=Yf(:,4)';

X_3Pts=affpt_to_pt3x3(affptx);  Y_3Pts=affpt_to_pt3x3(affpty);
% 将X和Y组合为一个矩阵用于输入：[Xpt1;Ypt1;Xpt2;Ypt2;Xpt3;Ypt3]
Xp = X_3Pts([1:2,4:5,7:8],:)';   
Yp = Y_3Pts([1:2,4:5,7:8],:)';  



