% -----------------------------------------
% Cluster permutation analysis - Within group results
% Anna Blasi Sept 2018
% Borja Blanco Jan 2022
% -----------------------------------------

clear; close all; clc

% Initialize variables for CLUSTER BASED PERMUTATION
% All 3-channel, triangle-shape combinations

LH = [1 2 3; 1 2 4; 1 3 4; 4 2 3; 4 2 6; 4 3 5; 4 5 6; 4 6 7; 4 5 7; ...
    7 6 5; 7 6 8; 7 5 9; 7 8 9; 7 8 10; 7 9 10; 10 8 9; 10 8 25; 10 9 11; ...
    10 11 25; 10 25 26; 10 11 26; 26 25 11];

RH = [12 13 14; 12 13 16; 12 14 16; 16 14 13; 16 13 18; 16 14 15; 16 18 15; ...
    16 18 17; 16 15 17; 17 18 15; 17 18 19; 17 15 20; 17 19 20; 17 19 22; ...
    17 20 22; 22 19 20; 22 19 24; 22 20 21; 22 24 21; 22 24 23; 22 21 23; ...
    23 21 24];

VEC_COMB = [LH; RH];

[NCOMB, sizeROI] = size(VEC_COMB);

% Initialize variables for storing data
obs_val_HbO(1:NCOMB)  = zeros;
P_HbO(1:NCOMB)  = zeros;

obs_val_HbR(1:NCOMB)  = zeros;
P_HbR(1:NCOMB)  = zeros;

NPerms = 1000;
nch=26;

% Load data, in format participants x channels
% Participant data represents the average value in the time window of interest
% for each participant, each channel and each condition
load('/Users/borjablanco/Library/CloudStorage/OneDrive-UniversityofCambridge/BB_documents/BASIS/STAARS_backup/STAARS/BB_analysis/analysis_072021/group_results/data4clusterPerm_tfce.mat')

% Select data/condition to analyze (within group)
ntl = 23; % typical likelihood participants
nel = 58; % elevated likelihood participants

% Run once for each group and condition (see commented example below)
% Visual social condition
%HbO = HbO_perm_s (1:ntl, :);
%HbR = HbR_perm_s (1:ntl, :);

% Vocal - non-vocal contrast
%HbO = HbO_perm_v(1:ntl, :) - HbO_perm_n(1:ntl, :);
%HbR = HbR_perm_v(1:ntl, :) - HbR_perm_n(1:ntl, :);

h = waitbar(0,'Computing Cluster Permutation');

for n = 1:NCOMB % number of channel combinations
    waitbar(n/NCOMB, h)

    % Select cluster for analysis
    current_ROI = VEC_COMB(n,:);
    
    [obs_val_HbO(n), P_HbO(n), obs_val_HbR(n), P_HbR(n)] = ClusterPermutationAnalysis_Basis(HbO, HbR, current_ROI, NPerms);

end
close(h)

% Store results
ClusterPermutationResults_HbO = table(VEC_COMB, obs_val_HbO', P_HbO','VariableNames',{'VEC_COMB' 'obs_val_HbO' 'P_HbO'});
ClusterPermutationResults_HbR = table(VEC_COMB, obs_val_HbR', P_HbR','VariableNames',{'VEC_COMB' 'obs_val_HbR' 'P_HbR'});

writetable(ClusterPermutationResults_HbO,['clusterPerm_10_16_EL_S_HbO' '_nperm1K']);
writetable(ClusterPermutationResults_HbR,['clusterPerm_10_16_EL_S_HbR' '_nperm1K']); 


