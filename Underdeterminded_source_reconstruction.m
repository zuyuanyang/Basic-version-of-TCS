 function  [ estS]                            =Underdeterminded_source_reconstruction(xw,estH,Nsources,m,stop)
%==========================================================================
K                                               =size(xw,3);
sT                                              =K/2+1;
P                                               =size(xw,2); 
T                                               =nchoosek(Nsources,m);
allcomb                                         =nchoosek(1:Nsources,m);
combinatedestS                            =zeros(Nsources,P,T,sT);
  p                                                =0.8;  %% Lp norm
%==========================================================================
       P                                        =size(xw,2);
range                                           =1:P;
estS                                            =zeros(Nsources,P,sT);
%==========================================================================
X                                               =permute(xw,[2,3,1]);
threshold                                       =10^3;
for k                                           =1:sT   
      hf                                        =squeeze(estH(:,:,k));
        if Nsources>m    
                 for t=1:T 
                 subHf                          =squeeze(estH(:,allcomb(t,:),k));
                if cond(subHf)<=threshold
                 index                          =allcomb(t,:);
                 combinatedestS(index,:,t,k)    =subHf\squeeze(xw(:,:,k)); 
                  elseif cond(subHf)>threshold
                  index                         =allcomb(t,:);    
                  combinatedestS(index,:,t,k)   =(subHf+10^-5.*randn(m,m))\squeeze(xw(:,:,k)); 
                  end
                 end
        end
end
%==========================================================================
 for ff                               =1:sT
     X                               =squeeze(xw(:,:,ff));
     hf                              =squeeze(estH(:,:,ff));
 %=========================================================================
 for jj                              =range 
          XW                         =squeeze(X(:,jj));
%==============================Initilization of sources============================================ 
for iter                            =1
  Lamda                             =rand(T,1);         
 Lamda                              =Lamda./sum(Lamda);%
 Start                              =zeros(Nsources,1);
 for kk                              =1:T
    weight                           =Lamda(kk);
   Start                             =Start+weight.*(combinatedestS(:,jj,kk,ff)); 
 end
%==========================================================================
[y,error]                           =complexaugmentedfocuss(XW,hf,m,Nsources,stop,p,Start);
end
%==========================================================================
 estS(:,jj,ff)                        =y;
 end
 end
estS                               =permute(estS,[2,3,1]);
estS(:,K:-1:sT+1,:)         =conj(estS(:,2:(sT-1),:));
 %%======================================================================== 












                 %=========================================================
%                      range                     =(d-1)*navg+1:d*navg;
%                      estS                      =squeeze(combinatedestS(index,range,t,k));
%                      cov                       =estS*estS';
%                    COV(index,index)            =cov;         
%                    error                       =squeeze(Rx(:,:,d,k))-squeeze(estH(:,:,k))*COV*squeeze(estH(:,:,k))';
%                    cov
%==========================================================================
%                    for i=1:m
%                      flag                      =allcomb(t,i);
%                      EstS(i,:)                 =squeeze(combinatedestS(flag,range,t,k));  
%                      location                  =index(i);
%                   estsita(location,d,t,k)      =(EstS(i,:)*EstS(i,:)');  
%                   end
%                   %%estsita(:,d,t,k)             =estsita(:,d,t,k);%./sum(estsita(:,d,t,k));
%                  combinationsita(t)            =max(estsita(:,d,t,k));
%                  end
%                    [val,loc]                   =max(combinationsita); 
%                  combinationsubspace(:,d,k)    =allcomb(loc,:);
%                  end          
%==========================================================================  
%         elseif Nsources<=m   %%% Overderdetermined case
%             hf                                 =estH(:,:,k);       
%            combinatedestS(:,:,k)               =pinv(hf)*squeeze(xw(:,:,k));
%         end     
% end    
%==========================================================================

