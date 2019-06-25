function     [eigenvectors,gamma,similarity]     =Fast_density_Clustering(percent,K,eigenvectors)
method                           =2; % 1. no scaling problem, 2, scaling calculation
%%=========================================================================
 for f                           =1:K/2+1 
%%========================================================================= 
blocknumber                       =size(eigenvectors{f},2);
%==========================================================================
referenceblock                    =1:blocknumber;
  for t                           =1:blocknumber
      index                       =(t-1)*blocknumber+1: t*blocknumber;
     y(index)                     =repmat(t,blocknumber,1);
  end  
   z                              =repmat(referenceblock',blocknumber,1);
%%=========================================================================
  [distance,~]                     =distance_calculation(eigenvectors{f}.',method);
   similarity{f}                   =distance;
     x                             =reshape(similarity{f},blocknumber^2,1); 
   xx(:,1)                         =y;
   xx(:,2)                         =z;
   xx(:,3)                         =x;
   [Gamma,Rho,Delta]               =originalcluster_dp(xx,percent);
 %---------------------------------------------------------------------------------
 Len                                             =size(Rho,2);
 %%========================================================================
  [~,Loc]                                       =find(Rho>=1); 
 %%======================================================================== 
  newGamma{f}                                   =Gamma(Loc);
  newRho                                        =Rho(Loc);
  newDelta                                      =Delta(Loc);
     eigenvectors{f}                                =eigenvectors{f}(:,Loc); 
   similarity{f}                                =similarity{f}(Loc,Loc);
 %=========================================================================  
   [~,newloc]                                    =sort(newGamma{f},'descend'); 
 %========================================================================= 
   Index{f}                                      =newloc;
    gamma{f}                                 =newGamma{f}(newloc);
   rho{f}                                        =newRho(newloc);
   delta{f}                                      =newDelta(newloc);
   eigenvectors{f}                          =eigenvectors{f}(:,newloc); 
   Similarity{f}                                =similarity{f}(newloc,newloc);
   xx=[];
   x=[];y=[];z=[];
 end
end