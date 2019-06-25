%--------------------------------------------------------------------------
% load speech/music sources
%--------------------------------------------------------------------------
if strcmp(channeltype, 'RIR')==1  
  [convolvedsource,mixture,origS,origH]      =RIRgeneration(sourcetype,RT60,Nsources,m,K);
end
%---------------------------- Music only ----------------------------------
 