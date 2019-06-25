function   [features ]      =Eigenvector_extraction(Rx,K)
%%=========================================================================
m                         =size(Rx,1);
 length                 =size(Rx,3); 
eigenvalue          =zeros(m,length,K/2+1);
 for f                    =1:K/2+1
    for b                 =1:length
        [U,E]             =eig(squeeze(Rx(:,:,b,f)));
        features{f}(:,b)   =U(:,end);
    end
   end
end
%%=========================================================================