% IEM_tutorial_EEG_advanced.m
%
% Now that you've spent some time learning the basics of how to run an IEM
% analysis on a small set of data, let's explore some more advanced aspects
% of the analysis you'll need to understand if you want to run an IEM-based
% experiment yourself.
%
% This tutorial will cover:
%
% - kinds of EEG data to use: focus on time course of alpha power
%
% - run-wise cross-validation: compute reconstructions for each trial by
%   looping over runs (leave-one-run-out). We'll also discuss other
%   cross-validation schemes. 
%
% - balancing trials within CV folds - due to artifact rejection, etc, you
%   sometimes end up with lopsided #'s of trials over feature space. I'll
%   show an example for how to address this, and we'll see how much it
%   impacts the data quality
%
% - Converting channel responses into 'reconstructions' - we'll weight the
%   channel sensitivity profile for each channel by its estimated
%   activation. This is a 'smoothed' version of the channel response
%   function (CRF), and can be 'exactly' coregistered across trials (not
%   just based on bin). This will also let us more fully explore how
%   different basis sets impact our analyses
%
% - Quantifying reconstructions based on their 'fidelity' (Sprague et al,
%   2016) and/or their 'slope' (Foster et al, 2016), the two most
%   commonly-used metrics right now. We'll also discuss when it's
%   appropriate to do curve-fitting to quantify 'amplitude' and 'size'
%
% - Computing reconstructions over time to uncover neural representations
%   at the ~ms timescale. We'll go over a few ways to do this, and discuss
%   how to compare results from each method. 
%
% - Understanding why it's important to compare reconstructions across
%   conditions (or timepoints) computed using a 'fixed' encoding model
%   (Sprague et al, 2018)
%
% - Evaluating the 'significance' of reconstructions: we'll generate a null
%   distribution of reconstructions against which we'll compare our actual
%   data. This is an important step to ensure any structure in the data is
%   not due to chance. We'll discuss these procedures at the group level as
%   well.


% where is the EEG data saved?
%eeg_root = '/Volumes/data/FosterEtAl2016/Data';
%subj_num = 10; % 3 extremely smooth; 8 is good, 10 good (delayed onset) - 10 is significant.. (750:1250)


% a few 'global' variables
pbins = 0:45:(360-45); % center of position bins
twindow = [-500 2000]; % the time window we'll consider for analyses

eeg_sampling_rate = 250; % recorded eeg data at 250 Hz


%% first - load data
%load(sprintf('%s/%02.f_WM_EEG.mat',eeg_root,subj_num));
addpath util/;
load data/EEG_advanced.mat;

% <document data>


% let's look at the 'raw' data - pick a channel (XXX) and sort trials by
% position bin and plot ERP; look at mean delay period potential vs
% position bin

elec_to_plot = 7;  % subj 8: elec_to_plot = 7
delay_window = [750 1250]; % for delay-period comparisons, which tpts?
delay_tpts = tpts >= delay_window(1) & tpts <=delay_window(2); % which tpts we'll use for 'delay' analyses

figure; subplot(1,2,1); hold on;
pu = unique(c_all(:,2)); % get all position indices
pos_colors = hsv(length(pu)); tmp_mu = nan(length(pu),1); tmp_std = nan(length(pu),1);
for pp = 1:length(pu)
    thisidx = c_all(:,2)==pu(pp) & ~excl_all; % which trials we average
    plot(tpts,squeeze(mean(dt_all(thisidx,elec_to_plot,:),1)),'Color',pos_colors(pp,:));
    tmp_mu(pp) = mean(mean(dt_all(thisidx,elec_to_plot,tpts>=delay_window(1)&tpts<=delay_window(2)),3),1);
    tmp_std(pp) = std(mean(dt_all(thisidx,elec_to_plot,tpts>=delay_window(1)&tpts<=delay_window(2)),3),[],1);
end
xlim(twindow);
xlabel('Time (ms)');
ylabel('Potential (\muV)');
title(sprintf('Electrode %s',chan_labels{elec_to_plot}));

subplot(1,2,2);
hold on; plot(pbins,tmp_mu,'k-');
for pp = 1:length(pu)
    plot(pbins(pp)*[1 1],tmp_mu(pp)+[-1 1]*tmp_std(pp),'-','Color',pos_colors(pp,:));
    plot(pbins(pp),tmp_mu(pp),'o','MarkerSize',10,'Color','w','MarkerFaceColor',pos_colors(pp,:))
end
xlabel('Position bin (center, \circ)');
ylabel('Mean delay potential (\muV)');
title(sprintf('%i to %i ms',delay_window(1),delay_window(2)));
xlim([-45 360]);

%% filter EEG data

% set up our filter properties

filt_band = [8 12]; % filter from 8-12 Hz

% NOTE: you could set this up to look at different freqs - left as exercise
% for reader...

% use EEGFILT to band-pass the data, then get instantaneous power via
% Hilbert transform

dbp_all = nan(size(dt_all)); % band-passed data
df_all = nan(size(dt_all));  % freq data
for cc = 1:size(dt_all,2)
    % band-pass data using eeg_filt
    dbp_all(:,cc,:) = eegfilt(squeeze(dt_all(:,cc,:)),eeg_sampling_rate,filt_band(1,1),filt_band(1,2));
    
    % hilbert-transform and compute power (abs(H(X)).^2)
    df_all(:,cc,:)  = abs(hilbert(eegfilt(squeeze(dt_all(:,cc,:)),eeg_sampling_rate,filt_band(1,1),filt_band(1,2)).').').^2; % instantaneous power calculated here for induced activity.
end

%% inspect filtered data!

% let's look at an example trial, channel:
tnum = 51; chnum = 5;
figure;
hold on;
% plot raw signal
plot(tpts,squeeze(dt_all(tnum,chnum,:)));
% plot band-pass signal
plot(tpts,squeeze(dbp_all(tnum,chnum,:)));
% plot abs(Hilbert-transform) - envelope
plot(tpts,squeeze(df_all(tnum,chnum,:)).^0.5);
% plot power from Hilbert
plot(tpts,squeeze(df_all(tnum,chnum,:)));

title(sprintf('Trial %i, Channel %i',tnum,chnum));
legend({'Raw','Filtered','Envelope','Power'});


% Similar to above, let's look at the same channel sorted by position
figure; subplot(1,2,1); hold on;
tmp_mu = nan(length(pu),1); tmp_std = nan(length(pu),1);
for pp = 1:length(pu)
    thisidx = c_all(:,2)==pu(pp) & ~excl_all; % which trials we average
    plot(tpts,squeeze(mean(df_all(thisidx,elec_to_plot,:),1)),'Color',pos_colors(pp,:));
    tmp_mu(pp) = mean(mean(df_all(thisidx,elec_to_plot,tpts>=delay_window(1)&tpts<=delay_window(2)),3),1);
    tmp_std(pp) = std(mean(df_all(thisidx,elec_to_plot,tpts>=delay_window(1)&tpts<=delay_window(2)),3),[],1);
end
xlim(twindow);
xlabel('Time (ms)');
ylabel(sprintf('%i to %i Hz power (\\muV^2)',filt_band(1),filt_band(2)));
title(sprintf('Electrode %s',chan_labels{elec_to_plot}));

subplot(1,2,2);
hold on; plot(pbins,tmp_mu,'k-');
for pp = 1:length(pu)
    plot(pbins(pp)*[1 1],tmp_mu(pp)+[-1 1]*tmp_std(pp),'-','Color',pos_colors(pp,:));
    plot(pbins(pp),tmp_mu(pp),'o','MarkerSize',10,'Color','w','MarkerFaceColor',pos_colors(pp,:))
end
xlabel('Position bin (center, \circ)');
ylabel('Mean delay power (\muV^2)');
title(sprintf('%i to %i ms',delay_window(1),delay_window(2)));
xlim([-45 360]);

% What can we deduce from comparing the alpha-filtered delay-period
% responses compared to mean delay-period responses (first figure; above)? 
% - try changing which electrodes, delay periods you plot


%% Build spatial IEM (same as in fMRI dataset!)
%
% Motivated by the spatial sensitivity of delay-period alpha responses in
% many electrodes, let's build a spatial IEM, which projects alpha activity
% from EEG electrode space into polar angle coordinates, so that we can
% directly average 'representations' across all trials. 
%
% This is overall the same procedure as fMRI, though SNR issues may require
% some adjustments to analysis. A few valid analysis choices involve:
% 
% - selecting 'selective' channels based on training data, or using
%   reported channels in the literature (e.g., Wolff et al, 2017)
% 
% - ensuring perfectly equated representation of each location bin in
%   training set: due to artifact rejection, it may be necessary to perform
%   model estimation using average exemplars within position bins rather
%   than the actual positions. Out of an abundance of caution, it would
%   also be necessary to ensure the same # of trials contribute to each
%   position bin in the training set, and shuffle this a few times to be
%   sure results don't come from accidental trial groupings
%
% - training on 'average' delay period alpha, then test on each timepoint?
%   (etc)



n_chan = 8; % # of channels, evenly spaced around the screen
chan_centers = linspace(360/n_chan,360,n_chan);

% evaluate basis set at these
angs = linspace(-179,180,360);

% TODO: replace!!!!
myb_orig = build_basis_polar_mat(angs,chan_centers);

figure; plot(angs,myb_orig,'LineWidth',1.5);
xlabel('Polar angle (\circ)');ylabel('Channel sensitivity'); title('Basis set (information channels)');
xlim([-180 180]);set(gca,'XTick',-180:90:180,'TickDir','out','Box','off');


% and, using our c_all variable, which contains the exact angle viewed on
% each trial, predict channel responses
%
% start by making a matrix n_trials x 360, with 1's and 0's, which we can
% multiply by our basis_set (b_orig) to generate predicted channel
% responses

% first, we'll need to wrap around angles greater than 180 to -179
c_all(c_all(:,1)>180,1) = c_all(c_all(:,1)>180,1)-360;

% sub2ind(size(stim_mask),1:size(stim_mask,1),find(c_all(:,1)==angs))
stim_mask = zeros(size(c_all,1),length(angs));
for tt = 1:size(c_all,1)
    % for now, we're going to assume the representation 'should be' perfect
    % - the feature mask we're projecting into channel space will be a
    % delta function. it's also possible to try a 'blob' here too.
    stim_mask(tt,angs==c_all(tt,1)) = 1;
    clear idx;
end

% HINT: here's a faster way to do this w/ newer versions of matlab:
% stim_mask = c_all(:,1)==angs;

% predict channel responses for each trial
% (when using a delta stimulus, this is equivalent to shifting a basis
% function and evaluting it at the channel centers - but sometimes your
% model may be more complicated than a delta stimulus, so it's good to be
% in the habit of writing it like this)
X_all = stim_mask * myb_orig;


% let's check out the predicted channel responses:
figure; 
% first, a single trial
whichtrial_C = 150; 
subplot(1,2,1);
hold on;  chan_colors = lines(n_chan);
plot(chan_centers,X_all(whichtrial_C,:),'k-','LineWidth',1.5);
for cc = 1:n_chan
    plot(chan_centers(cc),X_all(whichtrial_C,cc),'ko','MarkerSize',10,'MarkerFaceColor',chan_colors(cc,:))
end
plot([1 1]*mod(c_all(whichtrial_C),360),[0 1],'k--','LineWidth',2);
xlabel('Channel center (\circ)');
ylabel('Predicted channel response');
title(sprintf('Trial %i: %i\\circ',whichtrial_C,c_all(whichtrial_C,1)));
set(gca,'TickDir','out','XTick',0:90:360);


% all trials
subplot(1,2,2); hold on;
imagesc(chan_centers,1:size(c_all,1),X_all); axis ij;
plot([0 360],whichtrial_C+[-0.5 0.5],'r-');
xlim([22.5 360+22.5]);ylim([0.5 size(c_all,1)+0.5]);
xlabel('Channel center (\circ)');
ylabel('Trial');
title(sprintf('All trials (rank = %i)',rank(X_all)));
set(gca,'TickDir','out','XTick',0:90:360);


%% Train/test IEM (full delay period) - leave-one-run-out
% 
% With our encoding model and predicted channel responses in hand, now we
% can estimate channel weights and reconstruct WM representations from each
% channel. Since we only have one task here, we'll also perform
% leave-one-run-out cross-validation (to start with). We'll also explore a
% couple of other options for training/testing IEMs (including that used by
% Foster et al, 2016)

chan_resp = nan(size(X_all)); % fill this in with estimated channel responses

IEM_delay = [750 1250]; % time window to use for delay analyses (750:1250)
IEM_delay_tpts = tpts >= IEM_delay(1) & tpts <= IEM_delay(2);

% NOTE: here, you could change which data we're using! (electrodes, raw,
% etc)
delay_data = mean(df_all(:,:,IEM_delay_tpts),3);

% since we're doing leave-one-run-out (LORO), we'll use the run label to
% sort trials into training/testing sets
ru = unique(r_all);
for rr = 1:length(ru)
    trnidx = r_all~=ru(rr) & ~excl_all; % train using all 'included' trials except testing run
    tstidx = r_all==ru(rr); % for now, reconstruct with all trials (can drop excluded trials later)
    
    trndata = delay_data(trnidx,:);
    tstdata = delay_data(tstidx,:);
    trnX    = X_all(trnidx,:);
    
    w_hat = trnX \ trndata; % estimate channel weights
    
    chan_resp(tstidx,:) = tstdata/w_hat;
    
    clear w_hat trndata trnX tstdata trnidx tstidx;
end

%% Plot 'raw' reconstructed channel responses

% now we have channel responses computed on each trial
% note that these are 'raw' - different positions will result in different
% profiles of reconstructed channel response functions

% so let's average all trials of matching position bins and plot them, just
% to see if reconstructions track remembered position approximately
figure;
hold on;  pbin_plot = pbins; pbin_plot(pbin_plot==0)=360;
for pp = 1:length(pu)
    plot(chan_centers,mean(chan_resp(c_all(:,2)==pu(pp)&~excl_all,:),1),'-','LineWidth',2,'Color',pos_colors(pp,:));
    plot([1 1]*pbin_plot(pu(pp)),[0 1],'--','LineWidth',1.5,'Color',pos_colors(pp,:));
end
xlim([22.5 382.5]);
set(gca,'XTick',0:90:360,'TickDir','out');
xlabel('Channel center (\circ)');
ylabel('Reconstructed channel response (a.u.)');
title(sprintf('Delay-period reconstructions: %i to %i ms',IEM_delay(1),IEM_delay(2)));

%% Convert channel responses to 'reconstructions'
% 
% weight channel profile by channel response on each trial to get a smooth
% reconstruction, then shift & align all reconstructions together

% there are two ways we can do this (I'll show both)

% 1) circularly shift weighted channel profile to align all trials
recons_raw = chan_resp * myb_orig.';
recons_aligned = nan(size(recons_raw));

% we want to adjust so that each position is set to 0
for tt = 1:size(c_all,1)
    shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
    recons_aligned(tt,:) = circshift(recons_raw(tt,:),-1*shift_by);
end

% 2) adjust channel centers so weighted sum is aligned (this is more
% important for 2d models)
% TODO


% plot aligned reconstructions: all indiv trials and average
figure;
subplot(3,1,[1 2]);
imagesc(angs,1:sum(~excl_all),recons_aligned(~excl_all,:)); axis ij;
ylabel('Trial (after exclusion)');
xlabel('Angle (relative to target, \circ)');
title('Aligned delay-period reconstructions');
xlim([-180 180]);
set(gca,'TickDir','out','XTick',-180:90:180,'Box','off');

subplot(3,1,3); % mean across trials
hold on;
plot(angs,mean(recons_aligned(~excl_all,:),1),'k-','LineWidth',2);
% also plot +- SEM over trials (not best metric, but gives some sense of
% errorbars)
plot(angs,mean(recons_aligned(~excl_all,:),1) + [-1;1]*std(recons_aligned(~excl_all,:),[],1)/sqrt(sum(~excl_all)),'k--');
xlim([-180 180]);
set(gca,'XTick',-180:90:180,'TickDir','out');
xlabel('Angle (relative to target, \circ)');
ylabel('Reconstructed channel response (a.u.)');
title('Average over all trials');


%% how many trials are excluded from each bin? did this hurt us?

% now that we're pretty happy we have nice reconstructions during the delay
% period, we can/should check whether this could possibly be due to having
% more/less trials in a given bin on average than other bins

figure;
histogram(pbins(c_all(excl_all==0,2)),8);

% to me, this looks basically ok - but let's implement a version of the
% analysis that, within each position bin (on each CV fold) uses the same
% number of trials per bin


%% try IEM with fixed # of trials per bin
%
% this is almost identical to the leave-one-run-out CV above, except we're
% going to limit the training set a bit (so lots of this is copy/paste)

chan_resp_balanced = nan(size(X_all)); % fill this in with estimated channel responses

% NOTE: for simplicity, I'm going to just use the first n trials of each
% bin, where n is the number of trials in the bin with the fewest trials on
% that CV fold. Ideally, you'd repeat this procedure 10ish times, drawing a
% different set of trials each time. See Foster et al, 2016 paper & code
% for details on this procedure

for rr = 1:length(ru)
    trnidx = r_all~=ru(rr) & ~excl_all; % train using all 'included' trials except testing run
    tstidx = r_all==ru(rr); % for now, reconstruct with all trials (can drop excluded trials later)
    
    % look in c_all(:,2) at trnidx for each pos bin, see which has min
    % count
    incl_trials_bin = nan(length(pu),1); % how many trials we include for each bin
    binned_trials = cell(length(pu),1); binned_X = cell(length(pu),1);
    for pp = 1:length(pu)
        incl_trials_bin(pp) = sum(c_all(trnidx,2)==pu(pp));
        binned_trials{pp} = delay_data(trnidx & c_all(:,2)==pu(pp),:); 
        binned_X{pp} = X_all(trnidx & c_all(:,2)==pu(pp),:);
    end
    
    % within each bin, we'll only include this many
    incl_trials = min(incl_trials_bin); % note: can also do this w/ cellfun...
    binned_trials = cellfun(@(xx) xx(1:incl_trials,:), binned_trials,'UniformOutput',false);
    binned_X      = cellfun(@(xx) xx(1:incl_trials,:), binned_X,'UniformOutput',false);
    trndata = vertcat(binned_trials{:});
    trnX    = vertcat(binned_X{:});
    
    % loop over position bins, find trnidx trials within this bin, and add
    % the first incl_trials to trndata
    tstdata = delay_data(tstidx,:);
    
    w_hat = trnX \ trndata; % estimate channel weights
    
    chan_resp_balanced(tstidx,:) = tstdata/w_hat;
    
    clear w_hat trndata trnX tstdata trnidx tstidx;
end

% recenter, etc:
% 1) circularly shift weighted channel profile to align all trials
tmp_raw = chan_resp_balanced * myb_orig.';
recons_aligned_balanced = nan(size(tmp_raw));

% we want to adjust so that each position is set to 0
for tt = 1:size(c_all,1)
    shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
    recons_aligned_balanced(tt,:) = circshift(tmp_raw(tt,:),-1*shift_by);
    clear shift_by;
end
clear tmp_raw;

% now plot: we'll re-plot the original chan_resp's, and add the new
% balanced ones as dashed lines
% - TODO: try with balanced plotting set too

figure; subplot(3,1,[1 2]);
hold on;  pbin_plot = pbins; pbin_plot(pbin_plot==0)=360;
for pp = 1:length(pu)
    plot(chan_centers,mean(chan_resp(c_all(:,2)==pu(pp)&~excl_all,:),1),'-','LineWidth',2,'Color',pos_colors(pp,:));
    plot(chan_centers,mean(chan_resp_balanced(c_all(:,2)==pu(pp)&~excl_all,:),1),':','LineWidth',2,'Color',pos_colors(pp,:));
    plot([1 1]*pbin_plot(pu(pp)),[0 1],'--','LineWidth',1.5,'Color',pos_colors(pp,:));
end
xlim([22.5 382.5]);
set(gca,'XTick',0:90:360,'TickDir','out');
xlabel('Channel center (\circ)');
ylabel('Reconstructed channel response (a.u.)');
title(sprintf('Delay-period reconstructions: %i to %i ms',IEM_delay(1),IEM_delay(2)));

% and the average across all trials, aligned, smoothed
subplot(3,1,3); hold on;
plot(angs,mean(recons_aligned(~excl_all,:),1),'k-','LineWidth',2);
plot(angs,mean(recons_aligned_balanced(~excl_all,:),1),'k:','LineWidth',2);
xlim([-180 180]);
set(gca,'XTick',-180:90:180,'TickDir','out');
xlabel('Angle (relative to target, \circ)');
ylabel('Reconstructed channel response (a.u.)');
title('Average over all trials');
legend({'All trials','Balanced'});


%% Quantifying reconstructions: 'fidelity' (Sprague et al, 2016)
%
% Now that we have reconstructions, what do we do with them? Visualizing
% them as channel response functions and/or spatial reconstructions (as
% plotted above) is certainly useful - but it can be even more useful to
% distill the graded response profile into one or more parameters
% describing it. Here, we'll go over a few different ways to quantify these
% reconstructions (focusing in this tutorial on 1-d reconstructions, but
% similar principles apply to 2-d reconstructions, see Sprague & Serences,
% 2013).
%
% We'll start simple: let's compute a metric that is big when the
% reconstruction is, on average, peaked in the 'correct' position, and
% hovers around zero when there's no information. There are two ways to do
% this, both with substantial precedent in the literature. 
%
% First, let's compute what I call "FIDELITY" - this amounts to the
% 'strength' of the reconstruction in the 'direction' of the stimulus on
% that trial. Because, in this tutorial (and with color/orientation/polar
% angle stimuli), we're dealing with stimuli spanning a circular space, we
% can essentially just take the circular mean of the reconstruction. That
% is, if we plot the reconstruction in polar coordinates, the mean value of
% that plot can be considered a point-wise estimate of the reconstruction.
% The direction in which that mean vector points can be used to infer the
% represented feature value on that trial, and its length the 'strength'.
% But, of course the direction can be 'wrong' - so instead of reporting the
% 'strength' (vector length), we report the 'fidelity' - the vector length
% projected onto the 'correct' feature value. 

% Let's start with some examples. We'll plot reconstructions as above, averaged across
% each position bin, in polar coordinates. (Note that, for simplicity, I'm
% going to use the first set of reconstructions we computed - not the ones
% computed with 'balanced' training sets. Feel free to swap out data
% variables below to see how things change!)

% (note: we're using binned averages, not perfectly-aligned reconstructions
% here, for simplicity & visualization. you'd actually want to use the
% perfectly-aligned data, as we see below)

figure; ax1 = subplot(1,2,1,polaraxes); hold on;
for pp = 1:length(pu)
    thisidx = ~excl_all & c_all(:,2)==pu(pp);
    
    % plot the reconstruction in polar coords
    polarplot(deg2rad(angs),mean(recons_raw(thisidx,:),1),'-','Color',pos_colors(pp,:),'LineWidth',1.5);
    
    clear thisidx;
end
title('Aligned reconstructions (binned)');

% now let's focus on one bin (bin 2)
ax2 = subplot(1,2,2,polaraxes);  hold on; which_bin = 6;
title(sprintf('Example bin (%i)',which_bin));


% first, draw the same as above for 'background'

% get the trials we're using:
thisidx = ~excl_all & c_all(:,2)==which_bin;

polarplot(deg2rad(angs),mean(recons_raw(thisidx,:),1),'-','Color',pos_colors(which_bin,:),'LineWidth',1.5);

% and draw the mean of this position bin:
polarplot([1 1]*deg2rad(pbins(which_bin)),[0 1],':','Color',pos_colors(which_bin,:),'LineWidth',1);

% we want to do the vector average; let's add some example vectors
which_angs = (-135:45:180)-22; % (offsetting them from bin center)
% (note: these are JUST FOR VISUALIZATION) - we'll compute the real vector
% average below with all angles of the reconstruction. EXERCISE: is that
% really necessary? how few angles can we get away with? what are the
% constraints on which angles must be included?

for aa = 1:length(which_angs)
    % figure out which bin of recons_raw to average for each plotted angle
    angidx = find(angs==which_angs(aa)); % note: could replace w/ round, etc, for non-matching angles....
    polarplot([1 1]*deg2rad(which_angs(aa)),[0 1]*mean(recons_raw(thisidx,angidx),1),'-','LineWidth',2,'Color',[0.3 0.3 0.3]);
end

% now, we compute a vector sum of ALL points on the reconstruction
%
% but remember, this is a CIRCULAR variable - we can't just take the mean
% of recons_raw of these trials, so what do we do? 

% first, convert average reconstruction to X,Y
[tmpx,tmpy] = pol2cart(deg2rad(angs),mean(recons_raw(thisidx,:),1));

% now, take the mean of X,Y
mux = mean(tmpx); muy = mean(tmpy); clear tmpx tmpy;

% and convert back to polar coordinates:
[muth,mur] = cart2pol(mux,muy); % NOTE: muth is in RADIANS here

% and draw:
polarplot(muth*[1 1],mur*[0 1],'-','LineWidth',4,'Color',pos_colors(which_bin,:));

% Hopefully you can see that the average (thick colored line) aligns very
% nicely with the 'true' polar angle (thin dashed line) of this bin. So
% let's quantify this by projecting the vector mean (thick line) onto a
% unit vector in the direction of the bin center. Projection is just the
% dot product of the vector mean (mux,muy)<dot>(TrueX,TrueY), and TrueXY
% are just cosd(th) and sind(th), where theta is the bin center:
this_fidelity = dot([mux muy],  [cosd(pbins(which_bin)) sind(pbins(which_bin))]);
text(-pi/2,2,sprintf('Fidelity = %.03f',this_fidelity),'HorizontalAlignment','center');


% EXERCISE: does the order of operations here matter? can you compute
% fidelity on each trial then average, or do you need to compute fidelity
% on the average?



% just for fun, let's plot the vector averages of all trials color coded by
% their bin:
% all_mux = mean(recons_raw.*cosd(angs),2); all_muy = mean(recons_raw.*sind(angs),2);
% [all_muth,all_mur] = cart2pol(all_mux,all_muy);
% figure; polaraxes; hold on;
% for pp = 1:length(pu)
%     thisidx = ~excl_all & c_all(:,2)==pu(pp);
%     polarplot(all_muth(thisidx),all_mur(thisidx),'.','Color',pos_colors(pp,:));
% end

% So - we know how to project the vector average of the reconstruction onto
% the unit vector in the 'correct' direction, which lets us quantify
% fidelity. But, we've also computed 'aligned' reconstructions (indeed, we
% almost always do this in most analyses). In that case, the projection
% onto the 'correct' position amounts to projecting onto the x axis - so we
% can simplify our calculation of fidelity when using aligned data:
% F = mean( r(theta) * cos(theta) )  where r(theta) is 0-centered
% reconstruction value

% let's make a function that does this:
compute_fidelity = @(rec) mean( rec .* cosd(angs) ,2);

% compute fidelity for each trial
all_fidelity_delay = compute_fidelity(recons_aligned);

% plot distribution of fidelity across all trials
figure;
subplot(1,2,1);
histogram(all_fidelity_delay(~excl_all));
xlabel('Fidelity');
title('Fidelity across all trials');
set(gca,'TickDir','out');

% plot average fidelity (+std dev) for each position bin
subplot(1,2,2); hold on;
for pp = 1:length(pu)
    thisidx = ~excl_all & c_all(:,2)==pu(pp);
    thismu  = mean(all_fidelity_delay(thisidx));
    thiserr = std(all_fidelity_delay(thisidx))/sqrt(sum(thisidx));
    
    plot(pbins(pp)*[1 1],thismu + thiserr*[-1 1],'-','LineWidth',1.5,'Color',pos_colors(pp,:));
    plot(pbins(pp), thismu, 'wo','MarkerSize',8,'MarkerFaceColor',pos_colors(pp,:))
    clear thismu thiserr;
end
xlabel('Position bin center (\circ)');
ylabel('Mean fidelity');
title('Fidelity for each position bin');
tmpylim = get(gca,'YLim');
ylim([0 tmpylim(2)]);
xlim([-22.5 315+22.5]);
set(gca,'TickDir','out','XTick',0:90:360);

% what do you notice about the fidelity metric? Which positions are most
% poorly represented? Try using different position bins for the fidelity
% example computation above (change "which_bin")

%% Quantifying reconstructions: 'slope' (Foster et al, 2016)

% TODO


%% Quantifying reconstructions: curve fits (Ester et al, 2013)

% this typically works best for quantifying averages across trials, so
% we'll just do that here. For shuffling-type analyses like below, these
% fits can often be nonsensical, so it's also typically not the best way to
% establish data consistency/reliability. Instead, to get error bars of fit
% metrics, sometimes it can be good to bootstrap over trials. 



%% Reconstructions through time

% This is probably what you've all been waiting for! We're here because we
% think EEG is a cool tool, and often its most highly-regarded asset is its
% ability to measure neural activity at a very rapid timescale, much faster
% than can be afforded by techniques relying on indirect markers like
% BOLD fMRI. 
%
% There are a bunch of ways we can leverage this power. The most common,
% simplest, and most likely to work is to apply the IEM procedure at each
% timepoint during the trial in-turn. That is, to compute weights using a
% subset of trials at timepoint t, then reconstruct using activity at that
% same timepoint t. Let's give it a shot:

% which time range do we want to try reconstructing? 
recon_tpt_range = [-250 1750]; % from just before delay to end


% make a variable indexing into 3rd dim of data corresponding to
% interesting delay timepoints:
tptidx = find(tpts>=recon_tpt_range(1) & tpts<=recon_tpt_range(2));
% if we like, we can subsample:
ss_amt = 1; % how much we subsample - we take every n'th sample (1 = all samples, 2 = 1/2, 4 = 1/4, etc)
tptidx = tptidx(1:ss_amt:end);

% let's add an option to average neighboring timepoints/smooth (if we set
% this to the samplikng period, we won't be averaging at all)
avgwindow = 1000/eeg_sampling_rate; % average samples within this/2 on each side of tpt (ms)
%avgwindow = 50;

% make a variable to fill with reconstructions at each timepoint
chan_resp_t = nan(size(c_all,1),n_chan,length(tptidx)); % trial x channel x tpts

tpts_recon = tpts(tptidx); % when plotting reconstruciton timecourses, use these

for tt = 1:length(tptidx)
    
    % select which tpts to average for this iteration:
    % all timepoints within avgwindow/2 of tpts(tptidx(tt))
    this_tpt_idx = abs(tpts-tpts(tptidx(tt))) <= avgwindow/2; 
    
    % we can use this to index into our original data variable: df_all
    this_delay_data = mean(df_all(:,:,this_tpt_idx),3);
    
    for rr = 1:length(ru)
        trnidx = r_all~=ru(rr) & ~excl_all; % train using all 'included' trials except testing run
        tstidx = r_all==ru(rr); % for now, reconstruct with all trials (can drop excluded trials later)
        
        trndata = this_delay_data(trnidx,:);
        tstdata = this_delay_data(tstidx,:);
        trnX    = X_all(trnidx,:);
        
        w_hat = trnX \ trndata; % estimate channel weights
        
        chan_resp_t(tstidx,:,tt) = tstdata/w_hat;
        
        clear w_hat trndata trnX tstdata trnidx tstidx;
    end
    clear this_tpt_idx this_delay_data;
end

% recenter, etc:
tmp_raw = nan(size(chan_resp_t,1),length(angs),size(chan_resp_t,3));
for tt = 1:size(chan_resp_t,3)
    tmp_raw(:,:,tt) = chan_resp_t(:,:,tt) * myb_orig.';
end
recons_aligned_t = nan(size(tmp_raw));

% we want to adjust so that each position is set to 0
for tt = 1:size(c_all,1)
    shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
    recons_aligned_t(tt,:,:) = circshift(tmp_raw(tt,:,:),-1*shift_by);
    clear shift_by;
end
clear tmp_raw;

% plot the reconstructions as an image 
figure;
subplot(2,1,1); hold on;
imagesc(tpts_recon,angs,squeeze(mean(recons_aligned_t(~excl_all,:,:),1)));
%axis tight;
ylabel('Position (\circ)');
title(sprintf('Matched training; subsample: %ix, window: %i ms',ss_amt,avgwindow));
xlim(recon_tpt_range);ylim([angs(1) angs(end)]);
set(gca,'TickDir','out');

% and quantify: compute fidelity over time
fidelity_t = squeeze(compute_fidelity(recons_aligned_t)); % trials x tpts

subplot(2,1,2);
plot(tpts_recon,mean(fidelity_t(~excl_all,:),1),'k-','LineWidth',2);
xlabel('Time (ms)');
ylabel('Fidelity');
title('Fidelity through time (matched training)');
xlim(recon_tpt_range);
set(gca,'TickDir','out');


% and just because it's fun, a surface version
% TODO (from Josh's scripts)




%% Fixed encoding model across time
%
% Above, we recomputed the encoding model at each point in time, and this
% gave us nice, stable reconstructions through the entire delay. But how do
% we interpret that? If reconstructions are 'stronger' at one point in time
% compared to another, what does that mean? 
%
% It's hard to tell in this case. The reason is, the encoding model itself
% changes at each point in time, and so the space to whcih the IEM projects
% the electrode-space data is different! This tells us that there is
% information at each point in time, but does not say anything about the
% stability of this information - we cannot directly compare
% reconstructions at the beginning and end, except to say that there's
% more/less accessible information, not necessarily that the representation
% itself is waxing/waning. This is a subtle, but important point, and it
% applies to ALL comparisons of IEM results across conditions. If the model
% used to reconstruct across a comparison differs, we don't conclusively
% know whether the change in the reconstruction is due to a change in model
% fit, a change in the neural representation, or a combination of the two.
% This is a similar problem to seeing someone report a change in a ratio
% metric, without any report about the numerator or denominator
% individually - it is imposisble to know whether the numerator,
% denominator, or both changed across the manipulation of interest.
%
% One way around this is to use a 'fixed' encoding model - one estimated at
% a particular point in time, or using a separate set of data entirely
% (perhaps you have a localizer task, or a single-item set of trials, etc).
% Here, because there is only one 'condition', we'll estimate a fixed
% encoding model using the delay timepoints we were analyzing above. Then,
% we'll reconstruct using activity at each timepoint in turn. This will let
% us directly compare early- vs late-delay reconstructions, since they were
% computed using the same encoding model. 


% we'll use the same timepoints as above for reconstruction, but only a
% single epoch for training:
%trn_window = delay_window;
trn_window = IEM_delay;
trn_tpts = tpts >= trn_window(1) & tpts <= trn_window(2); delay_tpts; % from above - but experiment w/ different values!

chan_resp_t_fix = nan(size(c_all,1),n_chan,length(tptidx));


for rr = 1:length(ru)
    
    % on every cross-validation fold, train with a single set of exemplars
    % across trials averaged over trn_tpts

    trnidx = r_all~=ru(rr) & ~excl_all; % train using all 'included' trials except testing run
    tstidx = r_all==ru(rr); % for now, reconstruct with all trials (can drop excluded trials later)
    
    trndata = mean(df_all(trnidx,:,trn_tpts),3);
    trnX    = X_all(trnidx,:);
    
    % estimate IEM once per CV fold
    w_hat = trnX \ trndata; % estimate channel weights
    
    for tt = 1:length(tptidx)
        
        % select which tpts to average for this iteration:
        % all timepoints within avgwindow/2 of tpts(tptidx(tt))
        this_tpt_idx = abs(tpts-tpts(tptidx(tt))) <= avgwindow/2; 

        tstdata = mean(df_all(tstidx,:,this_tpt_idx),3);

        % use that estimated IEM (same across time!) to reconstruct each
        % timepoint
        chan_resp_t_fix(tstidx,:,tt) = tstdata/w_hat;
        
        clear this_tpt_idx tstdata;
    end
    clear w_hat trndata trnX trnidx tstidx;
end

% recenter, etc:
tmp_raw = nan(size(chan_resp_t_fix,1),length(angs),size(chan_resp_t_fix,3));
for tt = 1:size(chan_resp_t_fix,3)
    tmp_raw(:,:,tt) = chan_resp_t_fix(:,:,tt) * myb_orig.';
end
recons_aligned_t_fix = nan(size(tmp_raw));

% we want to adjust so that each position is set to 0
for tt = 1:size(c_all,1)
    shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
    recons_aligned_t_fix(tt,:,:) = circshift(tmp_raw(tt,:,:),-1*shift_by);
    clear shift_by;
end
clear tmp_raw;


% plot the reconstructions from fixed IEM as an image 
figure;
subplot(2,1,1); hold on;
imagesc(tpts_recon,angs,squeeze(mean(recons_aligned_t_fix(~excl_all,:,:),1)));
%axis tight;
ylabel('Position (\circ)');
title(sprintf('Fixed training (%i to %i ms)',trn_window(1),trn_window(2)));
xlim(recon_tpt_range);ylim([angs(1) angs(end)]);
set(gca,'TickDir','out');

% and quantify: compute fidelity over time
fidelity_t_fix = squeeze(compute_fidelity(recons_aligned_t_fix)); % trials x tpts

subplot(2,1,2); hold on;
plot(tpts_recon,mean(fidelity_t_fix(~excl_all,:),1),'k-','LineWidth',2);
plot(trn_window,[1 1]*-0.05,'k-','LineWidth',5);
xlabel('Time (ms)');
ylabel('Fidelity');
title('Fidelity through time (fixed training)');
xlim(recon_tpt_range);
set(gca,'TickDir','out');

% for completeness, let's compare directly:
figure;
hold on;
plot(tpts_recon,mean(fidelity_t(~excl_all,:),1),'-','LineWidth',2);
plot(tpts_recon,mean(fidelity_t_fix(~excl_all,:),1),'-','LineWidth',2);
plot(trn_window,[1 1]*-0.05,'k-','LineWidth',5); % training window
legend({'Matched training','Fixed training','Training window'});

ylabel('Fidelity');
xlabel('Time (ms)');
title('Comparing fidelity across IEM training regimes');


%% are these results by chance? shuffle trial labels at training
%
% to compute 'quality' of reconstructions on each iteration, we'll just use
% fidelity metric above (though slope, etc, also usable). Because we wrote
% functions to give us reconstruction metrics, we can just swap them in/out
% as we like! 
%
% How do we know the results we've shown above are 'real'? Especially for a
% single subject! A common way of doing this is removing structure from the
% 'training' dataset (typically by shuffling trial labels), computing an
% encoding model, and using that 'null' encoding model to reconstruct each
% trial of the test data (with their labels intact). Wash, rinse, repeat
% ~100-1000x. You can use this surrogate dataaset to compute any relevant
% stats, like fidelity/slope, and compare the 'real' value you see (at each
% timepoint even!) to this surrogate distribution to get a sense for the
% likelihood you'd see a real result by chance. 
%
% For simplicity, we'll just demonstrate this on the 'delay period' version
% of the analysis. But, it's quite easy to extend this approach to any of
% the temporal analyses we discussed (just computationally-intensive). 

niter = 500; % how many times do we shuffle?

% seed the random number generator
rng(2783829); % so we get the same result each time!

% initialize a few variables for 'null' reconstruction distributions
chan_resp_shuf = nan(size(c_all,1),n_chan,niter); % trials x channels x shuffling iterations
tmp_raw = nan(size(chan_resp_shuf,1),length(angs),niter);
for ii = 1:niter
    
    for rr = 1:length(ru)
        trnidx = r_all~=ru(rr) & ~excl_all; % train using all 'included' trials except testing run
        tstidx = r_all==ru(rr); % for now, reconstruct with all trials (can drop excluded trials later)
        
        trndata = delay_data(trnidx,:);
        tstdata = delay_data(tstidx,:);
        trnX    = X_all(trnidx,:);
        
        % shuffle rows of trnX before training
        trnX = trnX(randperm(size(trnX,1)),:);
        
        % now, there's no reliable correspondence between the predicted
        % response of each channel on a trial and the observed data on that
        % trial - any accidental structure that was there should be
        % removed!
        
        w_hat = trnX \ trndata; % estimate channel weights
        
        chan_resp_shuf(tstidx,:,ii) = tstdata/w_hat;
        
        clear w_hat trndata trnX tstdata trnidx tstidx;
    end

    tmp_raw(:,:,ii) = chan_resp_shuf(:,:,ii) * myb_orig.';
    
end

% now compute reconstructions, align


%tmp_raw = chan_resp * myb_orig.';
recons_aligned_shuf = nan(size(tmp_raw));

% we want to adjust so that each position is set to 0
for tt = 1:size(c_all,1)
    shift_by = c_all(tt,1); % if this is +, shift left (so use -1*shift_by)
    recons_aligned_shuf(tt,:,:) = circshift(tmp_raw(tt,:,:),-1*shift_by);
end
clear tmp_raw;

% and compute their fidelity
fidelity_shuf = squeeze(compute_fidelity(recons_aligned_shuf));

null_dist = mean(fidelity_shuf(~excl_all,:),1); % 100 'null' mean fidelity values
p = mean( null_dist >= mean(all_fidelity_delay(~excl_all)) );
figure;histogram(null_dist,10); hold on;
plot([1 1] * mean(all_fidelity_delay(~excl_all)), [0 15],'r--','LineWidth',2);
xlabel('Fidelity');
title(sprintf('Comparing true fidelity to shuffled (p = %0.03f)',p));



%% stable codes? how does info evolve over time in fixed; dynamic model

