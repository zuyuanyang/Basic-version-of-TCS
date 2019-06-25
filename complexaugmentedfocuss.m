function  [y,systemerror]         =complexaugmentedfocuss(X,A,m,n,stop,order,recover)
%==========================================================================
systemerror                         =zeros(1,stop);
for    i                                   =1:stop
       oldrecover                    =recover;
        diamat                        =(abs(recover)).^(2-order);
%===================================================================
      Diamat                          =diag(diamat);   
     tran                                =A*Diamat*A';
     if cond(tran)>10^2
      tran                               =A*Diamat*A'+10^-8*ones(m,m);   
     end
     recover                         =Diamat*A'*(tran\X);
       err                               =sum(abs(oldrecover-recover));
%         oldprob                   =abs(oldrecover).^order;
%         prob                      =abs(recover).^order;
%         err                       =sum(oldprob)-sum(prob); %
         systemerror(i)            =err;
  if  err<10^-4
    break;
  end
end
%=======================================================================
     y                           =recover;
%======================================================================== 