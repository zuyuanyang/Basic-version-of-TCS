clear; 
clc;
tic 
m                                    =2; % Number of sensors present
Nsources                         =3; % Number of sources present
K                                    =1024*2; %% FFT length 
sT                                   =K/2+1;
overlap                          =1/16;% K*ovlp must be an integer. Nonoverlap length.
%%=======Parameters useful for computation of the second-order statistics=====
 wind                            ='hanning'; 
 wlen                            =K;
 awin                            =hanning(wlen) ; % analysis window is a Hamming window
 navg                            =8;%       % Number of consecutive overlapping FFT frames that are used to compute the 
 %========================================================================= 
  sourcetype                       ='speech';
  channeltype                      ='RIR'; %%% RIR function
%=========================================================================
RT60                                       =0.1; %%% 100ms
%-------------------------RIR model----------------------------------------
cd build_mixtures
load_data;   
cd ..
close all;
signal_size                                              =size(mixture,1);
 Wmethod                                               ='whitening';
[xw,W]                                                   =FFTtransformation(mixture,overlap,K,m,Nsources,Wmethod,navg); %% STFT transformation
 [Rx,spectral]                                         =Corrf(xw,m,K,navg);%%% Calculate local covariance matrices
 [eigenvectors]                                    =Eigenvector_extraction(Rx,K); %%% Eigenvector extraction from local covariance matrices
 %%%=========================== Density based clustering =====================
 [est_H,eigenvectors]                            =Densitybased_clustering_method(K,eigenvectors,Nsources);
 if   strcmp(Wmethod,'whitening')==1
[xw,Rx,est_H,eigenvectors]                  =postwhitening(xw,sT,W,est_H,eigenvectors,Nsources,m,K,navg);
 end
  %%%=========================== K-weight clustering  adjustment=====================
[est_H]                                                =Weightclusteringmethod(eigenvectors,est_H,sT,m,Nsources);
%%%========= Sources reconstruction based on Lp norm measurement============
stop                                                     =50; %%%% Maximum iteration steps, e.g., stop =50; 
[estS]                                                  =Underdeterminded_source_reconstruction(xw,est_H,Nsources,m,stop);
%%%========= Permutation adjustment===================================
 [estsource,permuteS]                       =permutationprocess(Nsources,m,K,xw,est_H,estS); 
 %%%=========Recostruct sources are transformed into time domain===================================
separated                                          =reconstruct_output(estsource,K,overlap,signal_size,Nsources);
%%%========= Illustration===================================
resultdisplay(m,Nsources,mixture,origS,separated );