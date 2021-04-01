function [EV,REV,CorrM] = ExplainedVariance(Q_Pre, Q_Task, Q_Post)

% This script calculates explained variance (EV) and reverse EV
% the same manner it was calculated in Kudrimoti et al., 1999
%
% EV = ((R_task,post - R_task,pre*R_post,pre)/sqrt((1-R_task,pre.^2)(1-R_post,pre.^2))).^2;
%
% REV = ((R_task,pre - R_task,post*R_post,pre)/sqrt((1-R_task,post.^2)(1-R_post,pre.^2))).^2;
%
% It requires binned histograms of firing rates in full form (not sparse)
% of three intervals:
%  - activity during running epoch
%  - activity during preceding calm epoch
%  - activity during post-running calm epoch
%  
%  We calculate percentage of variance in post-running epoch that could be
%  explained by running epoch excluding the variance that comes from 
%  similarity between running epoch and pre-running epoch and
%  between pre- and post-running epochs
%  
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  
%     Inputs (not tsd!)
%     
%     - Q_Pre  -      binned histogram of firing rates
%                     during pre-running epoch (full! and better zscored)
%     - Q_Task  -     binned histogram of firing rates
%                     during running epoch (full! and better zscored)
%     - Q_Post  -     binned histogram of firing rates
%                     during post-running epoch (full! and better zscored)
%                     
%                     
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%      Outputs
%     
%     - EV  -      explained variance
%     - REV -      reverse expained variance
%                  (pre- and post-epochs are swopped)
%     - CorrM      Matrices of correlations for Pre-, Task and Post
%                  
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  
%  Coded by Dima Bryzgalov in MOBS team, Paris, France
%  06/03/2020

%% Check the form of the variables

% ToDo

%% Make sure the variables are not sparse
Q_Pre = full(Q_Pre);
Q_Task = full(Q_Task);
Q_Post = full(Q_Post);

%% Calculate the correlation matrices
CorrM.pre = corr(Q_Pre);
CorrM.task = corr(Q_Task);
CorrM.post = corr(Q_Post);

%% Find non-firing neurons in each epoch
idx_nonexist_pre = find(isnan(CorrM.pre(:,1)));
idx_nonexist_task = find(isnan(CorrM.task(:,1)));
idx_nonexist_post = find(isnan(CorrM.post(:,1)));

idx_toremove = unique(([idx_nonexist_pre; idx_nonexist_task; idx_nonexist_post])');

%% Remove non-firing neurons
CorrM.pre(idx_toremove,:) = [];
CorrM.pre(:,idx_toremove) = [];
CorrM.task(idx_toremove,:) = [];
CorrM.task(:,idx_toremove) = [];
CorrM.post(idx_toremove,:) = [];
CorrM.post(:,idx_toremove) = [];

%% Calculate correlation coefficients of correlation matrices

% R_task,post
R_TaskPost = corr2(CorrM.post,CorrM.task);
% R_task,pre
R_TaskPre = corr2(CorrM.task,CorrM.pre); 
% R_post,pre
R_PostPre = corr2(CorrM.post,CorrM.pre);

%% Calculate EV and REV

% Explained variance
EV = ((R_TaskPost-R_TaskPre*R_PostPre)/sqrt((1-R_TaskPre^2)*(1-R_PostPre^2)))^2;

% Reverse explained variance
REV = ((R_TaskPre-R_TaskPost*R_PostPre)/sqrt((1-R_TaskPost^2)*(1-R_PostPre^2)))^2;


end