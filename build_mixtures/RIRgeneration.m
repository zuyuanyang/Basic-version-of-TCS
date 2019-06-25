function  [convolvedsource,mixture,origS,origH]             =RIRgeneration(sourcetype,T60,I,J,numfreq)
%------------------------------------------------------------
if strcmp(sourcetype, 'speech')==1   
%========================================================================== 
  cd Sources
    [s1, Fs1]  =audioread('Iam_female_30s.wav');   
    [s2, Fs2]  =audioread('poem_male_30s.wav');    
    [s3, Fs3] = audioread('english_words_male_24s.wav'); 
    [s4, Fs4] = audioread('henry_theater_male_30s.wav');
   [s5, Fs5]  =audioread('sentence_female_28s.wav');
    [s6, Fs6] = audioread('numbers_male_30s.wav'); 
   cd ..
end
%==========================================================================
    A          =40;
 len          =Fs1*10;
%==========================================================================
       S1        =s1(1:len,1);
       S2        =s2(1:len,1);
       S3        =s3(1:len,1);
       S4        =s4(1:len,1); 
    S5        =s5(1:len,1);
  S6        =s6(1:len,1); 
%==========================================================================  
   origS(:, 1) = A*S1/norm(S1);       
   origS(:, 2) = A*S2/norm(S2);
   origS(:, 3) = A*S3/norm(S3); 
   origS(:, 4) = A*S4/norm(S4);  
  origS(:, 5) = A*S5/norm(S5); 
 origS(:, 6) = A*S6/norm(S6);   
     Fs          =Fs1;
%---------------------------------------------------------------------------------------------
% room dimensions (x,y,z)
rm=[5, 5, 2.3];
% reflection coefficients of the walls between -1 and 1 (e.g. -0.9 for large T60 and -0.1 for small T60)
% coefficient to control number of virtual sources 
V=rm(1)*rm(2)*rm(3);
S=rm(1)*rm(2)+rm(2)*rm(3)+rm(1)*rm(3);
S=2*S;
estalu=0.161*V/S/T60;
% ref=-sqrt(1-estalu)
ref=-(1-estalu);
%%=========================================================================
nvirt =40;
% Microphone coordinates (x,y,z)
mic1=[3 1 1.6]; mic2=[3 1.5  1.6]; mic3=[3 2 1.6]; mic4=[3 2.5 1.6]; %%%%%%%+(0.05*(test-1)) 1+(0.05*(monte-1)) 
mic5=[3 3 1.6]; mic6=[3 4 1.6]; mic7=[3 6 1.6]; mic8=[3 2 1.6];%%%    +(0.05*(test-1))
mic_mat=[mic1;mic2;mic3;mic4;mic5;mic6;mic7;mic8];
% source locations (x,y,z)
%%src1=[2.5+0.5*(monte-1) 1 1.6]; src2=[2.5+0.5*(monte-1) 1.4 1.6]; src3=[2 5 1.6]; src4=[2 3 1.6]; %*(test-1)
src1=[2   1 1.6]; src2=[2  1.4 1.6]; src3=[2 1.8 1.6]; src4=[2 2.2 1.6]; %%1+(0.4*(monte-1))
src5=[2 2.6 1.6]; src6=[2 4 1.6]; src7=[2 6 1.6]; src8=[2 2 1.6];
 src_mat=[src1;src2;src3;src4;src5;src6;src7;src8]; 
 
% Create impulse response for each (source-sensor) pair
% Create impulse response for each (source-sensor) pair
IR_all              =cell(J,I);
% H                  =newchannel.A;
% squeeze(H(1,1,1:2000))
for so              =1:I
    for m           =1:J
        h           = rir(Fs, mic_mat(m,:), nvirt, ref, rm, src_mat(so,:));
        IR_all{m,so}=h;
    end
end
%---------------------------------------------------------------------------
for i=1:J
for j=1:I
origH(i,j,:)              =fft(IR_all{i,j}(:),numfreq);% 
end
end
%%=========================================================================
%%%%visualize Room Impulse response between source 1 and microphone 1
figure;
Len=size(IR_all{1,1},2);
plot((1:Len)./Fs1,IR_all{1,1}); %(1/Fs)*
%axis([0 2000./Fs  min(IR_all{1,1}) max(IR_all{1,1})]);
axis([0 0.25 min(IR_all{1,1}) max(IR_all{1,1})])
%length(IR_all{1,1})
title('Impulse Response between source 1 microphone 1')
xlabel('Time (in seconds)');
ylabel('Magnitude');
for mic=1:J
   for so =1:I
       convolved(:, mic, so) = fftfilt(IR_all{mic, so}, origS(:, so));
    end
end
   mixture                          =sum(convolved(:,1:J,1:I),3); 
   for i=1:I
   convolvedsource(:,i)                         =squeeze(convolved(:,1,i)); 
   end