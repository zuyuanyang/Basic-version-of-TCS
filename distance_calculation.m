function      [distance,block]   =distance_calculation(X,method)
block                     =size(X,1);
%--------------------------------------------------------------------------
% for i                     =1:block
%     if real(X(i,1))<0
%     X(i,:)                =-X(i,:);
%     end
%     X(i,:)                =X(i,:)./norm(X(i,:));           
% end
%--------------------------------------------------------------------------
percent                =5;
%==========================================================================
pp                      =0.9; 
if method               ==1
  Inner                 =(X*X');
   distance             =sqrt(ones(block,block)-real(Inner));%-conj(Inner));
elseif method           ==2
      Inner             =(X*X');
 distance               =(ones(block,block)-Inner.*conj(Inner));%-conj(Inner)); 
 %=========================================================================
 ii                      =1;
 for i                   =1:block
     for j               =i+1:block
 xx(ii)                  =distance(i,j);
 ii                      =ii+1;
     end
 end
 %=========================================================================
%  position                =round(block^2*percent/100);
% sda                      =sort(xx);
% dc                       =0.1;%sda(position);
% xy                       =origIsomap(distance,2,[]);

 %distance                =1-exp(-distance.^pp);
 % distance                =exp(-distance./dc);
 
 % max(distance)
end
end