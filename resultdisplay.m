function           resultdisplay(m,n,mixture,S,estS)
%==========================================================================
figure
T                       =size(mixture(:,1),1);
for i=1:m
amplitude(i)            =max(abs(mixture(:,i)));
subplot(m,1,i)
plot(mixture(:,i));
 axis([0,T,-amplitude(i), amplitude(i)]);
  xlabel('Time (in sec.)');
   ylabel('Magnitude');
end
  xlabel('Observed mixture signals');
  
  
Or1                      =S(:,1);
Or2                      =S(:,2);
estx1                    =estS(:,1);
estx2                    =estS(:,2);
T                        =size(S(:,1),1);
for i=1:n
amplitude(i)            =max(abs(S(:,i)));
end
%==========================================================================
figure
for i=1:n
subplot(n,1,i)
plot(S(:,i))
 axis([0,T,-amplitude(i), amplitude(i) ]);
  xlabel('Time (in sec.)');
   ylabel('Magnitude');
end
  xlabel('Original sources');
T                    =size(estS(:,1),1);
%==========================================================================
figure
for i=1:n
amplitude(i)            =max(abs(estS(:,i)));
subplot(n,1,i)
plot(estS(:,i))
 axis([0,T,-amplitude(i), amplitude(i)]);
   xlabel('Time (in sec.)');
   ylabel('Magnitude');
end
  xlabel('Estimated sources');
