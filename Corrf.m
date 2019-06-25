function    [Rx,spectral]                 =Corrf(xw,J,K,Navg)
%--------------------------------------------------------------------------
    sT                                  =K/2+1;
   P                                    =size(xw,2);
method                                  =1;
%--------------------------------------------------------------------------
 D                                     =Navg;
 T                                     =floor(P/D);
%-------------------------------------------------------------------------- 
if method ==1
origRx                                  =zeros(J^2*T,sT);
xw                                      =xw(:,:,1:sT);	%  JxP *sT
%%=========================================================================
%%Short-time de-mean processing
% for t=1:T
%         index                           =(t-1)*D+1:t*D;
%         for ff=1:sT
%            Calmean                      =mean(squeeze(xw(:,index,ff)),2);
%            xw(:,index,ff)               =xw(:,index,ff)-repmat(Calmean,1,D);
%         end
% end
%%=========================================================================
 % Compute the cross-correlation for each time-frequency point
Indc                                    = reshape((1:J)'*ones(1,J),J^2,1);		% [1:J .. 1:J]'
Indr                                    = reshape(ones(J,1)*(1:J),J^2,1);		% [1..1 2..2 J..J]'
Innerproduct                            = xw(Indc,:,:).*conj(xw(Indr,:,:));     % sT x J^2 x T
for t=1:T
        index                           =(t-1)*D+1:t*D;
       origRx(J^2*(t-1)+1:J^2*t,:)      =squeeze(sum(Innerproduct(:,index,:),2));
end
%--------------------------------------------------------------------------
Rx                                     =reshape(origRx,J,J,T,sT) ;
spectral                               =zeros(T,sT);
 for f =1:sT
     Rx_mat                            =reshape(Rx(:,:,:,f),J*J,T);
     for d=1:T
         compRx                        =abs(Rx_mat(:,d));
     spectral(d,f)                     =sum(compRx);
     end
%%-------------------------------------------------------------------------    
%      for t=1:T
%          [VV,DD]                      = eig(Rx(:,:,t,f));
%         sigma2_est_vec                = min(diag(DD));
%         Rx(:,:,t,f)                   =Rx(:,:,t,f)-min(sigma2_est_vec)*eye(J); % denoise
%      end   
%%-------------------------------------------------------------------------- 
 end  
 end  
end
%--------------------------------------------------------------------------
 