function [gamma,rho,delta]=originalcluster_dp(xx,percent)

ND=max(xx(:,2));
NL=max(xx(:,1));
if (NL>ND)
  ND=NL;
end
N=size(xx,1);

dist  =zeros(ND,ND);
for i=1:N
  ii=xx(i,1);
  jj=xx(i,2);
  dist(ii,jj)=xx(i,3);
  dist(jj,ii)=xx(i,3);
end
%%%%percent=8;
%%%%%fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);

position=round(N*percent/100);
sda=sort(xx(:,3));
dc=sda(position);

rho  =zeros(1,ND);

for i=1:ND
     if i>1 && i<ND
     index     =[1:i-1,i+1:ND];
     elseif i==ND
     index     =1:i-1;  
     elseif i==1
      index    =i+1:ND;   
     end
    % pp=1.8;
     pp=2;
     rho(i)   =sum(exp(-((dist(i,index).^pp)./(dc.^2))));    
end
%%=========================================================================
maxd                =max(max(dist));
[rho_sorted,ordrho]  =sort(rho,'descend');
%%=========================================================================
delta(ordrho(1))=-1;
pp =5;
nneigh(ordrho(1))=0;
for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
       delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
       % delta(ordrho(ii))= exp(-(dc^2/dist(ordrho(ii),ordrho(jj))^pp));
        nneigh(ordrho(ii))=ordrho(jj);
     end
   end
end
delta(ordrho(1))=1;%max(delta(:));
maxrho          =max(rho);
minrho          =min(rho);
maxdelta        =max(delta);
mindelta        =min(delta);
truedelta       =maxdelta-mindelta;
%===========================
%  rhotao         =minrho+0.1*(maxrho-minrho);
%  rho            =rho.*exp(-rhotao./rho);
% deltatao        =mindelta+0.4*truedelta;
%  delta          =delta.*exp(-(deltatao)./(delta));
%===========================
gamma              =rho.*delta;%
end