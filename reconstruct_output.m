function normestS =reconstruct_output(U,K,ovlp,signal_size,Nsources)
%%=========================================================================
    sizeU          =size(U);
    separated      =zeros(signal_size,sizeU(3)); 
 for p=1:sizeU(3)
         for n=1:K*ovlp:signal_size -K   
            separated(n:n+K-1,p)=separated(n:n+K-1,p)+real(ifft(squeeze(U(:,(n+K*ovlp-1)/(K*ovlp),p))));
         end
 end
    A                                                   =40;
   len                                                  =size(separated,1);
for i=1:Nsources
     Sep                                                =separated(1:len,i);
   normestS(:,i)                                     = A*Sep/norm(Sep);       
   Sep=[];
end