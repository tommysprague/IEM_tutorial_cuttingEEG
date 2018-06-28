% IEM_tutorial_basics.m
%
% First, let's go over a first set of basic concepts underlying the IEM
% method without all the overhead of processing EEG datasets. I've already
% set up data to be simple to use, and this is a high-fidelity fMRI dataset
% that has quite robust single-trial decoding performance. For our
% purposes, the dataset is basically identical to the EEG datasets we'll be
% working with later, but this one shoudl be a bit less cumbersome to deal
% with when 
%
% In this tutorial, we'll learn:
% 1) how to build a channel-based linear encoding model
% 2) how to use that encoding model, together with a behavioral task, to
%    predict how modeled channels should respond
% 3) how to compute the contribution of each channel to each measurement
%    (here, activity from EEG electrode)
% 4) how to use that estimated ENCODING MODEL to solve for channel
%    responses on a separate set of data
% 5) how to think about these reconstructions


%% Load data

% I've already minimally-processed this dataset so that we have a measured
% scalp activity pattern (alpha) on each trial. We'll cover ways to get
% this pattern later on. You can also try loading fMRI_basic.mat - the data
% is stored exactly the same way, and the datasets are interchangeable for
% our analyses here. The point is that any signal can be a good candidate
% for this analysis, so long as some of the assumptions (discussed below)
% more or less hold. 

load('data/EEG_basic.mat'); % or fMRI_basic.mat, if interested - the data is stored hte same, just different source!

% in your workspace, you should have:
%
% - data_all: one activity pattern on each trial from each measurement
%  (electrode or voxel) - n_trials x n_measurements
%
% - c_all: condition label for each trial. first column is the polar angle
%   remembered on that trial, in degrees (Cartesian coords, 45 = upper
%   right); second column is the position bin (1-8; beginning at 0 deg and
%   moving CCW, so bin 2 is centered at 45, 3 at 90, etc); n_trials x 2
%
% - excl_all: logical vector, contains 1 for each trial that was marked by
%   Foster et al (2016) for exclusion based on artifacts, etc; n_trials x 1
%
% - r_all: label for each trial as to which run it came from; n_trials x 1
%
% - chan_labels (for EEG): cell array of strings labeling each electrode;
%   n_electrodes+2 x 1 cell array (we don't include EOG channels in data
%   here)

% this is all we need! we'll go over a few suggestions for how to process
% EEG data for best use with IEM methods a bit later. 


%% Build encoding model

% The central assumption in the IEM analysis framework is that we can model
% the activity of an aggregate neural measurement (EEG electrode; fMRI
% voxel) as a linear combination of a discrete set of modeled 'information
% channels', which we sometimes call 'basis functions'. For now, let's
% think of those information channels as a set of neural populations - we
% know, based on decades of vision science, that neural populations in the
% visual system are TUNED to particular features. 
%
% What are some examples of these?
% - orientation
% - motion
% - spatial position <---- we'll focus on this one, since it's easiest
% 
% A very important part of the IEM approach is to build SIMPLIFIED models -
% we know that spatial RFs in visual cortex tile the whole visual field,
% from fovea to periphery, etc. But do we need to model all known aspects
% of these neural populations? Especially for something like EEG, NO. We
% need to model the *relevant dimensions* for our task. In the datasets
% we're working with, stimuli were always presented along a ring at equal
% eccentricity - so we can only measure changes in neural responses as a
% function of polar angle - there is no point in trying to model different
% responses for different stimulus eccentrificites, since we never measured
% that! (see Sprague & Serences, 2013; 2014; 2016 for examples of this
% analysis across eccentricity as well)
%
% The point is: THE MODEL SHOULD REFLECT YOUR UNDERSTANDING OF NEURAL
% SENSITIVITY PROFILES AND THE CONSTRAINTS OF YOUR EXPERIMENTAL DESIGN
%
% At an extreme level: consider modeling sensitivity to color (e.g.,
% Brouwer & Heeger, 2009) for the datasets we're using here. The stimuli we
% used never spanned that dimension, so such modeling is futile! We'll see
% below some of the mathematical/theoretical constraints on these models. 
%
% Ok - back to modeling. 

% The IEM framework we're using in this tutorial is LINEAR: that means, for
% each electrode, we're making the direct assumption that
%
% B = W1 * C1 + W2 * C2 + ... + Wn * Cn
% 
% Where B is measured activity on each trial, Ci is the activity of the i'th
% modeled information channel on that trial, and Wi is the WEIGHT
% describing how strongly channel Ci contributes to the electrode in
% question. 





%% Use encoding model to predict channel responses





%% Use predicted channel responses to fit encoding model

% Requirements: at least as many unique positions as modeled channels


% First let's look at one electrode

% How well does the model fit? Compare predicted response to measured


% And we can do this for all electrodes in a single step, using linear
% regression (this step is univariate!)


% How do changes in encoding model properties (above) change the model
% fits? Try adjusting # of channels, channel width, etc. 



%% Use estimated encoding model to reconstruct channel responsees
