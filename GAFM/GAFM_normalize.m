

function  [X, Y, normal] =GAFM_normalize(x, y)

if size(x,2)==6

    n=size(x,1);
    m=size(y,1);

    %% 
    normal.xd=mean([x(:,1:2); x(:,3:4); x(:,5:6)]);   normal.yd=mean([y(:,1:2); y(:,3:4); y(:,5:6)]); %
    
    x=x-repmat(normal.xd,n,3);
    y=y-repmat(normal.yd,m,3);

    %% 
    normal.xscale=sqrt(sum(sum(x.^2,2))/3/n); 
    normal.yscale=sqrt(sum(sum(y.^2,2))/3/m);

    X=x/normal.xscale;
    Y=y/normal.yscale;
else
    disp('error')
end
end

