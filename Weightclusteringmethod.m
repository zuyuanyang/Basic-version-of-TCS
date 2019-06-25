function    [est_H]                     =Weightclusteringmethod(features,est_H,sT,m,Nsources) 
iter_max                                =6;
 percent                                =3;
%==========================================================================   
for f                                   =1:sT   
for iter                                =1:iter_max 
       T                                 =size(features{f},2);   
 for i=1:Nsources
      X                                 =squeeze(est_H(:,i,f))';
      Y                                 =features{f};
      Inner                             =(X*Y);
    xx                                  =(ones(1,T)-Inner.*conj(Inner));  
 position                                =round(T*percent/100);
sda                                      =sort(xx);
dc                                       =sda(position);   
 weight                                  =exp(-(xx./dc).^2);
Weight                                   =repmat(weight,m,1);
newfeatures                             =Weight.*features{f}; 
%%%%=======================================================================
Feature                                =[real(newfeatures);imag(newfeatures)];
cov                                    =Feature*Feature';    
[U,E]                                  =eig(cov);
est_H(:,i,f)                           =U(1:m,end)+U(m+1:2*m,end)*1i;
%%%%%%=====================================================================
 weight                                =[];
 newfeatures                           =[];
 Feature                               =[];
%%%%=======================================================================
est_H(:,i,f)                           =est_H(:,i,f)./norm(est_H(:,i,f));
end
end
end
%=================================
for f                                    =1:sT 
     for i                              =2:m
     est_H(i,:,f)                    =est_H(i,:,f)./est_H(1,:,f);      
     end
     est_H(1,:,f)                  =ones(1,Nsources);
end
end
%==========================================================================
