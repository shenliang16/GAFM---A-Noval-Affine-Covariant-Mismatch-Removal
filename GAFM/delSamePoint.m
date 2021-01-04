function [loc,des,loc4]=delSamePoint(loc,des,loc4)
if isempty(loc)
   loc = []; des=[]; loc4=[];
   return; 
end

temp=tril(abs([loc(:,1)]-[loc(:,1)]')<0.25,-1)&tril(abs([loc(:,2)]-[loc(:,2)]')<0.25,-1);
loc(any(temp,2),:)=[];       
if nargin>1
    des(any(temp,2),:)=[];
else
    des = []; 
end
if nargin>2
    loc4(any(temp,2),:)=[];
else
    loc4 = [];
end

% n=size(loc,1);
% for i=1:n
%     n=size(loc,1);
%     if i<n
%         temp=abs(loc(i,1:2)-loc(i+1:min(i+10,n),1:2))<0.01;
%         index=all(temp,2);
%         index_temp=find(index>0);
%         if ~isempty(index_temp)
%             loc(index_temp+i,:)=[]; des(index_temp+i,:)=[];
%         end
%     end
% end

