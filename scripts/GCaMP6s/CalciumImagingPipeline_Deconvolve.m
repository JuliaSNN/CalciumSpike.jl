function [FDeconv,PropNeg] = CalciumImagingPipeline_Deconvolve(F,Fneu,EpisodeNum,FramesPerEpisode,dt,subtract_Neuropil,tau)

%%
% Neuropil subtraction
% DeltaF/F calculation
% Deconvolve
% Gaussian filtering
%
% INPUT
%   F: fluorescence from suite2P
%   Fneu : fluorescence of surrounding neuropil from suite2P
%   EpisodeNum : number of episodes in the recording, we estimate the
%   baseline episode by espisode
%   FramesPerEpisdode : number of frame sin each episode, we estimate the
%   baseline episode by espisode
%   dt : time bin duration
%   subtract_Neuropil: (boolean, mandatory argument, default is yes (1))
%   tau : tau for deconvolution (optional), default is 2
% OUTPUT
%  FDeconv : deoncolved, filtered and generally beautified deltaF/F
%  PropNeg : for each neuron, the f0 during each episode, this will allow
%  to estimate how often this value is <0. If its <0 that's bad news

%% Set the basic parameters
NeuropilCoeff = 0.7; %0.7 is from Chen et al, 2013, calibration on GCAMP6s
FiltWindow = 100; % the size of the window for filtering data before baseline estimation
f0Prctle = 20;% take lowest f0Prctle% of values as a more robust estimate of baseline, I think is was 10 for gif data

%% Neuropil subtration
if subtract_Neuropil
    disp('removing neuropil')
    FSub = F - NeuropilCoeff*Fneu;
else
    FSub = F;
end

if ~exist('tau')
   tau = 2;
end

%% DeltaF / F calculation
Start = 1;
DeltaFOverF = [];
PropNeg = [];

% each neuron is processed in parallel by applying the matlab functions
% filt, prectile etc to the matrix

disp('DeltaF / F calculation')
for rec = 1:length(FramesPerEpisode) % if there are different recordis with different nmuber fo frames
    for ep = 1:EpisodeNum(rec)
        ep
        % get the data from the episde
        temp                   = FSub(:,Start : Start+FramesPerEpisode(rec)-1)';
        
        % smooth the activity before extracing baseline
        baseline               = fn_filt(temp,FiltWindow,'lk',1);
        
        % extract the baseline value, instead of taking the minimum take
        % lowest 20% percentile to avoid crazy outliers
        fo                     = prctile(baseline,f0Prctle);
        
        % Get delatF/F
        DeltaFOverF  = [DeltaFOverF,(temp ./ repmat(fo,[size(temp,1) 1])  - 1)'];
       
        PropNeg = [PropNeg;fo];
        
        Start = Start+FramesPerEpisode(rec);
    end
end

PropNeg = nanmean(PropNeg<0,1);
%% Deconvolution and filtering single trials
% Deconvolution equation : dV/dt + V.dt/tau = firing
disp('Deconvolution')

for CellNumber = 1:size(DeltaFOverF,1) % # of cells
    temp  = squeeze(DeltaFOverF(CellNumber,:));
    % deconvolve
    temp = [diff(temp,1,2) 0]+[temp(:,1:end-1)  0]*dt/tau;
    % this line is from original JB code, but should probably be [diff(temp,1,2) 0]*tau/dt+[temp(:,1:end-1)  0];

    % filter (gaussianfilter.m: reconstructed in this same folder — the
    % original was lost, see gaussianfilter.m's docstring)
     FDeconv(CellNumber,:) = gaussianfilter(temp,1,1.5,1)'; % parameters from Deneux et al, 2017

end

end