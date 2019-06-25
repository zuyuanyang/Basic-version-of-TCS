function       [xw,W]    =FFTtransformation(x,overlap,T,M,N,Wmethod,navg)
%--------------------------------------------------------------------------
 [N,J]                           = size(x);
 for mic=1:J
     x(:,mic)                  =x(:,mic)-mean(x(:,mic));
 end
%%=======================================================================
[signal_size Nmixtures]  =size(x);
win                      =window(@hann,T);
%converting mixture into time frequency domain
X                       =zeros(T,ceil((signal_size-T)/(T*overlap)),Nmixtures);        %initialise buffer
for m=1:Nmixtures
    for k=1:T*overlap:signal_size - T      
        xw(:,(k+T*overlap-1)/(T*overlap),m) = fft(x(k:k+T-1,m).*win);
    end
end
 xw                     =permute(xw,[3,2,1]); %% M*Q*F
 %=========================================================================
 Q                     =floor(size(xw,2)/navg);
 for f                 = 1:T/2+1
       for b           =1:Q
       cov             =squeeze(xw(:,(b-1)*navg+1:b*navg,f))*squeeze(xw(:,(b-1)*navg+1:b*navg,f))'; 
       end
 end
 %================================================
  if   strcmp(Wmethod,'whitening')==1 
for f= 1:T
    cov         =(squeeze(xw(:,:,f))* squeeze(xw(:,:,f))') ;
   R_bar        = cov;
[U S]           =svd(0.5*(R_bar+R_bar'));
B               =U(:,1:M)*sqrt(S(1:M,1:M));
W(:,:,f)        =B;
xw(:,:,f)       =B\xw(:,:,f);
end
  elseif   strcmp(Wmethod,'whitening')==0
      R_whiten=[];
      W=[];
  end 
end
 
 