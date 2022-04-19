function [EV,REV] = ExplainedVariance_Behavior(mapPre, mapTask, mapPost)

% This function calculates explained variance (EV) and reverse EV
% of occupation maps. Inpired by reactivation metrics suggested in
% Kudrimoti et al., 1990
%
% EV = ((R_task,post - R_task,pre*R_post,pre)/sqrt((1-R_task,pre.^2)(1-R_post,pre.^2))).^2;
%
% REV = ((R_task,pre - R_task,post*R_post,pre)/sqrt((1-R_task,post.^2)(1-R_post,pre.^2))).^2;
%
% It requires occupation maps of equal size from:
%  - activity during epoch of interest (EOI)
%  - activity during period preceding EOI
%  - activity during post-EOI epoch
%  
%  We calculate percentage of variance in post epoch that could be
%  explained by EOI excluding the variance that comes from 
%  similarity between EOI and pre-epoch and
%  between pre- and post-epochs
%  
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  
%     Inputs
%     
%     - mapPre  -     occupation map during PreEpoch
%     - mapTask  -    occupation map during epoch of interest
%     - mapPost  -    occupation map during PostEpoch
%                     
%                     
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%      Outputs
%     
%     - EV  -      explained variance
%     - REV -      reverse expained variance
%                  (pre- and post-epochs are swopped)
%                  
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  
%  Coded by Dima Bryzgalov in MOBS team, Paris, France
%  Modified from ExplainedVariance.m
%  28/06/2021


%% Calculate correlation coefficients of correlation matrices

% R_task,post
R_TaskPost = corr2(mapPost,mapTask);
% R_task,pre
R_TaskPre = corr2(mapTask,mapPre);
% R_post,pre
R_PostPre = corr2(mapPost,mapPre);

%% Calculate EV and REV

% Explained variance
EV = ((R_TaskPost-R_TaskPre*R_PostPre)/sqrt((1-R_TaskPre^2)*(1-R_PostPre^2)))^2;

% Reverse explained variance
REV = ((R_TaskPre-R_TaskPost*R_PostPre)/sqrt((1-R_TaskPost^2)*(1-R_PostPre^2)))^2;

end