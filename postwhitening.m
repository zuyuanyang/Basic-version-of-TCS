function [xw,Rx,est_H,features]= postwhitening(xw,sT,W,est_H,features,Nsources,m,K,navg)
for f =1:sT
    features{f}                                             =W(:,:,f)*features{f}; 
  est_H(:,:,f)                                              =W(:,:,f)*est_H(:,:,f); 
    xw(:,:,f)                                               =W(:,:,f)*xw(:,:,f);
   for i=1:Nsources
       est_H(:,i,f)                                         =est_H(:,i,f)./norm(est_H(:,i,f));
         if  real( est_H(1,i,f)  )>=0
       est_H(:,i,f)   = est_H(:,i,f)  ;
         elseif real(est_H(1,i,f) )<0
       est_H(:,i,f)   =-est_H(:,i,f)  ;   
        end  
   end
   
   for j=1:size(features{f},2)
          features{f}(:,j)                                      =features{f}(:,j)./norm(features{f}(:,j)); 
                 if  real( features{f}(1,j) )>=0
        features{f}(:,j)   = features{f}(:,j) ;
                 elseif real( features{f}(1,j) )<0
        features{f}(:,j)  =- features{f}(:,j) ;   
                  end     
   end
  end
  [Rx,origRx]                                               =Corrf(xw,m,K,navg);
