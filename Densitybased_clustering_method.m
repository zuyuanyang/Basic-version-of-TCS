function    [estH,eigenvectors]                                    =Densitybased_clustering_method(F,eigenvectors,Nsources)
 percent                                                                      =6;
 sT                                                                              =F/2+1  ;
[eigenvectors,gamma,similarity]                                =Fast_density_Clustering(percent,F,eigenvectors);
 [est_clusternumber]                                                  =Clusters_detection(similarity,gamma,F) ;
  [estH ]                                                                      =Sort_eigenvectors(similarity,eigenvectors,est_clusternumber,F,Nsources);
for f                                                                             =1:F/2+1                  
[val,loc]                                                                       =sort(gamma{f},'descend');
eigenvectors{f}                                                            =eigenvectors{f}(:,loc);
end
end
%============================================================================
 function [estH ]                                                           =Sort_eigenvectors(similarity,eigenvectors,est_clusternumber,K ,estsourcenumber)
method                                                                        =1;
for f                                                                             =1:K/2+1                  
  binestsourcenumber                                                =est_clusternumber(f);
    index                                                                      =1:binestsourcenumber;
    buindex                               =binestsourcenumber+1:estsourcenumber;
     if method                             ==1
       estH(:,:,f)                         =eigenvectors{f}(:,1:estsourcenumber);   
     elseif method                         ==2
   P                                      =size(eigenvectors{f},2);
ratio                                     =1;
P1                                        =floor(P*ratio);
    if binestsourcenumber >=estsourcenumber
       distance                            =similarity{f}(1:binestsourcenumber,1:binestsourcenumber);
       for i=1:binestsourcenumber
           distance(i,i)=10^10;
       end
loc                                    =0;
range                                  =[];
newdistance(1)                           =10^10;
  for k                                 =2:estsourcenumber  
       [ newdistance(k),~]             =min(distance(k,1:binestsourcenumber));
  end
  [val,loc]                             =sort(newdistance,'descend');
   estH(:,:,f)                         =eigenvectors{f}(:,loc(1:estsourcenumber));  
    elseif binestsourcenumber <estsourcenumber
estH(:,1:binestsourcenumber,f)             =eigenvectors{f}(:,1:binestsourcenumber);  
               distance                    =similarity{f}(1:P1,1:P1);
loc                                    =0;
range                                  =1:binestsourcenumber;
  for k                                 =binestsourcenumber+1:estsourcenumber 
         newvolume                      =zeros(1,P1);   
       for i=k:P1
            if i~=range
           [ newvolume(i),~]            =min(distance(i,range));
            elseif i==range
           newvolume(i)                 =0;     
            end
       end
  [val,loc]                             =max(newvolume);
    estH(:,k,f)                         =eigenvectors{f}(:,loc); 
  range                                 =[range,loc];                       
  end
  buji                                  =setdiff(1:P1,range);
    end
     end
end
 end
%==========================================================
 function      [est_clusternumber]             =Clusters_detection(similarity,gamma,F)  
%--------------------------------------------------------------------------
sT                                             =F/2+1;
 est_clusternumber                                 =zeros(1,sT);
 cumparameter                                  =zeros(1,sT);
 est_ratio                                     =zeros(1,6);
%==========================================================================
for f                                                 =1:sT
       [K,value,sorte]                         =Gapbasedmethod(gamma{f});
       est_clusternumber(f)              =K; 
       SORTE{f}                               =sorte;
%========================================================================== 
      index                                          =1:K; 
      maximum                                   =floor(0.25*F);%%0.2*F   
     cumparameter(f)                         =(sum(sum(similarity{f}(index,index))))./K ;   
end
  [~,loc]                                            =sort(cumparameter,'descend');
 confidenceloc                                =loc(1:maximum);
%%=========================================================================
for i=1:6
     [~,loc]                                       =find(est_clusternumber(confidenceloc)==i);     
     est_ratio(i)                                =size(loc,2);
end
[~,ratioloc]                                   =max(est_ratio);
 end
%=================================
 function [K,newlada,sorte]  =Gapbasedmethod(newlada)
K                     =size(newlada,2);
range                  =[2:6];
sorte                  =zeros(1,6);
 for i                 =1:K-1  
    f(i)              =newlada(i)- newlada(i+1);
end
for j                =1:(K-1)
    landa(j)         =(1/(K-j))*sum(f(j:K-1));
end
for j                =1:(K-1)
    delta(j)          =(f(j)-landa(j))^2;
end
for j                =1:(K-1)
   sita(j)           =(1/(K-j))*sum(delta(j:K-1));
end
%%=========================================================================
for j                =1:(K-2)
    if sita(j)>0
  sorte(j)           =sita(j+1)/sita(j);
    elseif sita(j)==0
        sorte(j)     =1;  
    end
end
 [~,sorloc]          =find(sorte>0);
[val,loc]           =min(sorte(sorloc));
if val>0
outputv              =loc;
end
range                  =[2:5];
[sourceval,sourceloc]      =min(sorte(range));
K         =sourceloc+1;
 end


 function   [est_H,Adjustloc]           =Adjustchannelestimation(similarity,features,gamma,rho,sT,Nsources,bin_estnumber)
 %=========================================================================
for f                                  =1:sT
 if bin_estnumber(f)             ==Nsources
 est_H(:,:,f)                          =features{f}(:,1:Nsources);
 Adjustloc(:,f)                        =1:Nsources;
 %=========================================================================
 elseif bin_estnumber(f)               >Nsources  
 est_H(:,1,f)                           =features{f}(:,1);
  for i                                 =2:Nsources  
    index                               =[1:i]; 
    Sim                                 =similarity{f};
   [val,~]                              =min(Sim(i,1:i-1)); 
    gamma{f}(i)                         =val*rho{f}(i); 
  end
 [~,loc]                               =sort(gamma{f},'descend');
 est_H(:,2:Nsources,f)                 =features{f}(:,loc(2:Nsources));
 Adjustloc(:,f)                        =[1,loc(2:Nsources)]; 
%%==========================================================================
%    sourcenumber                        =bin_estnumber(f);  
%   index                               =1;
%   Wholeset                           =1:bin_estnumber(f); 
% for i                                 =2:Nsources  
%     gamma{f}                          =zeros(1,bin_estnumber(f)); 
%   Sim                                 =similarity{f}(index,:);  
%    Testindex                          =setdiff(Wholeset,index);
%   for t                              =1:size(Testindex,2)
%      [val,loc]                       =min(Sim(:,Testindex(t))); 
%      detla{f}(Testindex(t))          =val;
%      gamma{f}(Testindex(t))          =val*rho{f}(Testindex(t));
%   end
%  [val,loc]                          =max(gamma{f});  
%  index                              =[index,loc];
%  Sim                                =[];
% end
%  est_H(:,:,f)                         =features{f}(:,index);
%  Adjustloc(:,f)                       =index;    
%========================================================================== 
 elseif bin_estnumber(f)              <Nsources  
%   sourcenumber                        =bin_estnumber(f);  
%   est_H(:,1:sourcenumber,f)           =features{f}(:,1:sourcenumber);
%       index                           =1:sourcenumber; 
%   for i                                 =sourcenumber+1:size(similarity{f},2)   
%     Sim                                 =similarity{f};
%    [val,~]                              =min(Sim(index,i)); 
%     gamma{f}(i)                         =val*rho{f}(i); 
%   end
%  [~,loc]                               =sort(gamma{f},'descend');
%  est_H(:,sourcenumber+1:Nsources,f)    =features{f}(:,loc(sourcenumber+1:Nsources));
%  Adjustloc(:,f)                        =[1:sourcenumber,loc(sourcenumber+1:Nsources)]; 
%==========================================================================  
sourcenumber                        =bin_estnumber(f);    
  index                               =1:sourcenumber;
for i                                 =sourcenumber+1:Nsources   
    Wholeset                          =1:size(similarity{f},2);
    gamma{f}                          =zeros(1,size(similarity{f},2)); 
  Sim                                 =similarity{f}(index,:);  
   Testindex                          =setdiff(Wholeset,index);
  for t                              =1:size(Testindex,2)
     [val,loc]                       =min(Sim(:,Testindex(t))); 
     detla{f}(Testindex(t))          =val;
     gamma{f}(Testindex(t))          =val*rho{f}(Testindex(t));
  end
 [val,loc]                          =max(gamma{f});  
 index                              =[index,loc];
 Sim                                =[];
end
 est_H(:,:,f)                        =features{f}(:,index);
 Adjustloc(:,f)                      =index;
%==========================================================================
end
end
 end
 