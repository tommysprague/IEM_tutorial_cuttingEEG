function tcs_saveEEGdata(subs)
%
% save out EEG data like MRI data (trialwise)
%
% save out like fMRI data: trials x channels x timepoints
%
% want:
% - dt_all
% - c_all (polar angle of each position and bin) - pos angle seems to go positive
%   as CW from 0 deg cartesian, so we'll convert to same coord system as
%   MRI (+ang is CCW)
% - r_all (based on 64-trial blocks)
% - resp_all (response polar angle)
% - x_all (excluded trials)
%
% NOTE: posBin cetnered on 0, 45, 90, etc - where 45 is bottom right
%
% (adapted from SpatialEM.m from JJF and DWS, by TCS 6/25/2018)

subs = [1,2,3,5,7,8,9,10,12,13,14,15,16,18,20];
nSubs = length(subs);



% setup directories
root = pwd; out = 'AnalysisScripts';
%dRoot = [root(1:end-length(out)),'/Data/'];
%eRoot = [root(1:end-length(out)),'/EEG/'];
%bRoot = [root(1:end-length(out)),'/Behavior/'];

dRoot = fullfile(root,'Data');
eRoot = fullfile(root,'EEG');
bRoot = fullfile(root,'Behavior');

% parameters to set
%em.nChans = 8; % # of channels
%em.nBins = em.nChans; % # of stimulus bins
%em.nIter = 10; % # of iterations
%em.nBlocks = 10; % # of blocks for cross-validation
em.frequencies = [8 12]; % frequency bands to analyze
em.bands = {'Alpha'};
em.time = -500:4:2000; % time points of interest
em.window = 4;
em.Fs = 250;
em.nElectrodes = 20;

name = sprintf('WM_EEG.mat'); % name of files to be saved

% for brevity in analysis
% nChans = em.nChans;
% nBins = em.nBins;
% nIter = em.nIter;
% nBlocks = em.nBlocks;
freqs = em.frequencies;
times = em.time;
nFreqs = size(em.frequencies,1);
nElectrodes = em.nElectrodes;
nSamps = length(em.time);
Fs = em.Fs;

% Specify basis set
% em.sinPower = 7;
% em.x = linspace(0, 2*pi-2*pi/nBins, nBins);
% em.cCenters = linspace(0, 2*pi-2*pi/nChans, nChans);
% em.cCenters = rad2deg(em.cCenters);
% pred = sin(0.5*em.x).^em.sinPower; % hypothetical channel responses
% pred = wshift('1D',pred,5); % shift the initial basis function
% basisSet = nan(nChans,nBins);
% for c = 1:nChans;
%     basisSet(c,:) = wshift('1D',pred,-c); % generate circularly shifted basis functions
% end
% em.basisSet = basisSet; % save basis set to data structure

% Loop through participants
for s = 1:nSubs
    sn = subs(s);
    fprintf('Subject:\t%d\n',sn)
    
    % Grab data------------------------------------------------------------
    
    % Get position bin index from behavior file
    fName = [dRoot, '/', num2str(sn), '_MixModel_wBias.mat']; load(fName);
    em.posBin = beh.trial.posBin'; % add to fm structure so it's saved
    posBin = em.posBin;
    
    % Get EEG data
    fName = [eRoot, '/', num2str(sn), '_EEG.mat']; load(fName);
    eegs = eeg.data(:,1:20,:); % get scalp EEG (drop EOG electrodes)
    artInd = eeg.arf.artIndCleaned.'; % grab artifact rejection index
    %tois = ismember(eeg.preTime:4:eeg.postTime,em.time); nTimes = length(tois); % index time points for analysis.
    nTimes = size(eeg.data,3);
    
    excl_all = artInd; % 1 for exclude, 0 for keep
    c_all = [mod(-1*beh.trial.pos.',360) beh.trial.posBin.'];
    
    % posBin goes 1,2,3 from centered at 0, 45, 90, where +ang is CW from
    % cartesian 0 (right) - we want to flip it so that bin 2 is topright,
    % bin 3 is up, bin 4 is top left, etc. 
    c_all(c_all(:,2)~=1,2) = -1*c_all(c_all(:,2)~=1,2)+10; 
    
    
    resp_all = -1*beh.trial.err; % response error on each trial
    
    
    % create a list of run numbers
    % 64 trials per run, 15 runs
    r_all = reshape( repmat( (1:15),64,1 ),960,1 );
    
    tpts = eeg.preTime:(1000/eeg.sampRate):eeg.postTime;
    chan_labels = eeg.chanLabels;
    
    % Remove rejected trials
    %eegs = eegs(~artInd,:,:);
    %posBin = posBin(~artInd);
    
    em.nTrials = length(posBin); nTrials = em.nTrials; % # of good trials
    
    %----------------------------------------------------------------------
    
    % Preallocate Matrices
    %tf_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nChans); tf_total = tf_evoked;
    %C2_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nBins,nChans); C2_total = C2_evoked;
    %em.blocks = nan(nTrials,nIter);  % create em.block to save block assignments
    
    % Loop through each frequency
    %for f = 1:nFreqs
    %    tic % start timing frequency loop
        %fprintf('Frequency %d out of %d\n', f, nFreqs)
        
        % Filter Data
        %fdata_evoked = nan(nTrials,nElectrodes,nTimes);
        dt_all = nan(nTrials,nElectrodes,nTimes);
        ef_all = nan(nTrials,nElectrodes,nTimes); % filtered EEG
        
        dt_all = eegs;
        
     %   for c = 1:nElectrodes
            %fdata_evoked(:,c,:) =    hilbert(eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(f,1),freqs(f,2))')';
            
%            dt_all(:,c,:) = abs(hilbert(eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(f,1),freqs(f,2))')').^2; % instantaneous power calculated here for induced activity.
%            ef_all(:,c,:) = eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(f,1),freqs(f,2));
      %  end
    %end
    
    
    
    

        % Loop through each iteration
%         for iter = 1:nIter
%             
%             %--------------------------------------------------------------------------
%             % Assign trials to blocks (such that trials per position are
%             % equated within blocks)
%             %--------------------------------------------------------------------------
%             
%             % preallocate arrays
%             blocks = nan(size(posBin));
%             shuffBlocks = nan(size(posBin));
%             
%             % count number of trials within each position bin
%             clear binCnt
%             for bin = 1:nBins
%                 binCnt(bin) = sum(posBin == bin); 
%             end
%                       
%             minCnt = min(binCnt); % # of trials for position bin with fewest trials
%             nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
% 
%             % shuffle trials
%             shuffInd = randperm(nTrials)'; % create shuffle index
%             shuffBin = posBin(shuffInd); % shuffle trial order
%             
%             % take the 1st nPerBin x nBlocks trials for each position bin.
%             for bin = 1:nBins;   
%                 idx = find(shuffBin == bin); % get index for trials belonging to the current bin
%                 idx = idx(1:nPerBin*nBlocks); % drop excess trials
%                 x = repmat(1:nBlocks',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
%             end
%             
%             % unshuffle block assignment
%             blocks(shuffInd) = shuffBlocks;
%             
%             % save block assignment
%             em.blocks(:,iter) = blocks; % block assignment
%             em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
%             
%             %-------------------------------------------------------------------------
%             
%             % Average data for each position bin across blocks   
%             posBins = 1:nBins;
%             blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
%             blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
%             labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
%             blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
%             c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
%             bCnt = 1;
%             for ii = 1:nBins
%                 for iii = 1:nBlocks
%                     blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
%                     blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
%                     labels(bCnt) = ii;
%                     blockNum(bCnt) = iii;
%                     c(bCnt,:) = basisSet(ii,:);
%                     bCnt = bCnt+1;
%                 end
%             end
%      
%            for t = 1:nSamps
%                 
%                 % grab data for timepoint t
%                 toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
%                 de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
%                 dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
%                 
%                 % Do forward model
%                                 
%                 for i=1:nBlocks % loop through blocks, holding each out as the test set
%                     
%                     trnl = labels(blockNum~=i); % training labels
%                     tstl = labels(blockNum==i); % test labels
%                     
%                     %-----------------------------------------------------%
%                     % Analysis on Evoked Power                            %
%                     %-----------------------------------------------------%
%                     B1 = de(blockNum~=i,:);    % training data
%                     B2 = de(blockNum==i,:);    % test data
%                     C1 = c(blockNum~=i,:);     % predicted channel outputs for training data                 
%                     W = C1\B1;          % estimate weight matrix
%                     C2 = (W'\B2')';     % estimate channel responses
%                     
%                     C2_evoked(f,iter,t,i,:,:) = C2; % save the unshifted channel responses
%                     
%                     % shift eegs to common center
%                     n2shift = ceil(size(C2,2)/2);
%                     for ii=1:size(C2,1)
%                         [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                         C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                     end
%                     
%                     tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                     
%                     %-----------------------------------------------------%
%                     % Analysis on Total Power                             %
%                     %-----------------------------------------------------%
%                     B1 = dt(blockNum~=i,:);    % training data
%                     B2 = dt(blockNum==i,:);    % test data
%                     C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                     W = C1\B1;          % estimate weight matrix
%                     C2 = (W'\B2')';     % estimate channel responses
%                     
%                     C2_total(f,iter,t,i,:,:) = C2;
%                     
%                     % shift eegs to common center
%                     n2shift = ceil(size(C2,2)/2);
%                     for ii=1:size(C2,1)
%                         [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                         C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                     end
%                     
%                     tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                     %-----------------------------------------------------%
%                     
%                 end
%             end
%         end
%         toc % stop timing the frequency loop
%     end
    
    %fName = [dRoot, '/' , num2str(sn),name];
    fName = sprintf('%s/%02.f_%s',dRoot,sn,name);
    %em.C2.evoked = C2_evoked;
    %em.C2.total = C2_total;
    %em.tfs.evoked = tf_evoked;
    %em.tfs.total = tf_total;
    %em.nBlocks = nBlocks;
    save(fName,'dt_all','c_all','r_all','excl_all','tpts','chan_labels','-v7.3');
end
