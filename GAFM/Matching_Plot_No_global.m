
function Matching_Plot_No_global(Im11,Im22,loc1,loc2,Corresp)

% 统一大小
[m1,n1,D1]=size(Im11);  [m2,n2,D2]=size(Im22);  D=min(D1,D2); sizeMin=max([m1,n1],[m2,n2]);  
Im1=zeros([sizeMin,D]); 
Im2=zeros([sizeMin,D]);
if D1>D; Im11=rgb2gray(Im11); end;     Im1(1:m1,1:n1,:)=Im11; Im1=uint8(Im1);
if D2>D; Im22=rgb2gray(Im22); end;     Im2(1:m2,1:n2,:)=Im22; Im2=uint8(Im2);


interval = 20; WhiteInterval = uint8(255*ones(size(Im1,1), interval, D));
% figure('position',[120 120 1350 900]);
figure;
imagesc(cat(2, Im1, WhiteInterval, Im2)); 
colormap('gray'); hold on; cols1 = size(Im1,2);


line([loc1(Corresp(:,1),2)'; loc2(Corresp(:,2),2)'+cols1+interval], [loc1(Corresp(:,1),1)' ;  loc2(Corresp(:,2),1)'],'linewidth', 1, 'color','b' ) ;%[0,0.5,0.8]
axis equal; axis off;
end