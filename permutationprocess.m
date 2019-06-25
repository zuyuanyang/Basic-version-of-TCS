function      [estsource,perU]          =permutationprocess(Nsources,m,K,xw,estH,eS) 
 eS                                                  =permute(eS,[2,1,3]);
sT                                      =K/2+1; 
D                                       =size(xw,2);
for f                                   =1:K/2+1
    for i                               =1:Nsources
    Shat(:,f,i,1)                       =eS(f,:,i);
    end
    for jj=2:m
        for i=1:Nsources
            Shat(:,f,i,jj)              =eS(f,:,i)*estH(jj,i,f);
        end
    end
end
%==========================================================================    
estS                                  =squeeze(Shat(:,:,:,1));
%==========================================================================
  P                                    =size(estS,2);
sT                                     =K/2+1;
test                                   =1;
range                               =1:P;
%==========================================================================    
estS(:,K:-1:sT+1,:)                             =conj(estS(:,2:(sT-1),:));
[S_hat_p,perm_p]                                =perm_align(estS);
%==========================================================================
P                                               =size(S_hat_p,1);
estsource                                       =zeros(K,P,Nsources); 
U                                               =S_hat_p;
perU                                            =permute(U,[2,1,3]);
for i=1:Nsources
estsource(1:sT,:,i)                             =perU(1:sT,:,i);
estsource(K:-1:sT+1,:,i)                        =conj(perU(2:(sT-1),:,i));
end
%==========================================================================
%-------------------------------------------------------------------------