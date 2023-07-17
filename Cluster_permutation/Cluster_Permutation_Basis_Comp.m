% -----------------------------------------
% Cluster permutation analysis - Between group comparisons
% Anna Blasi Sept 2018
% Borja Blanco Jan 2022
% -----------------------------------------

clear; close all; clc

% Load data, in format participants x channels
% Participant data represents the average value in the time window of interest
% for each participant, each channel and each condition
load('/Users/borjablanco/Library/CloudStorage/OneDrive-UniversityofCambridge/BB_documents/BASIS/STAARS_backup/STAARS/BB_analysis/analysis_072021/group_results/data4clusterPerm.mat')

% Select data/condition to analyze (see commented example below)
% Visual social condition
%HbO = HbO_perm_s;
%HbR = HbR_perm_s;
 
nTL = 23;
nEL = 58;

% All 3-channel, triangle-shape combinations on BASIS layout
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
NPerms = 1000;
obs_val_HbO(1:NCOMB)  = zeros;
perm_vals_HbO = zeros(NPerms, NCOMB);

obs_val_HbR(1:NCOMB)  = zeros;
perm_vals_HbR = zeros(NPerms, NCOMB);

h = waitbar(0,'Computing Cluster Permutation');

for n = 1:NCOMB
    waitbar(n/NCOMB, h)

    % Select cluster for analysis
    current_ROI = VEC_COMB(n,:);
    
    [obs_val_HbO(n), perm_vals_HbO(:,n), obs_val_HbR(:,n), perm_vals_HbR(:,n)] = ClusterPermutationAnalysis_Basis_Comp(HbO, HbR, current_ROI, NPerms, nTL, nEL);

end
close(h)

%BB - p-val:
Tdist_HbO = sort(abs(perm_vals_HbO), 'ascend');

P_HbO(1:NCOMB) = zeros;
for i=1:NCOMB
    P_HbO(i) = (length(find(Tdist_HbO(:,i)>=abs(obs_val_HbO(i)))))/NPerms;
end

Tdist_HbR = sort(abs(perm_vals_HbR), 'ascend');

P_HbR(1:NCOMB) = zeros;
for i=1:NCOMB
    P_HbR(i) = (length(find(Tdist_HbR(:,i)>=abs(obs_val_HbR(i)))))/NPerms;
end

% Store results
ClusterPermutationResults_HbO = table(VEC_COMB, obs_val_HbO', P_HbO','VariableNames',{'VEC_COMB' 'obs_val_HbO' 'P_HbO'});
ClusterPermutationResults_HbR = table(VEC_COMB, obs_val_HbR', P_HbR','VariableNames',{'VEC_COMB' 'obs_val_HbR' 'P_HbR'});

writetable(ClusterPermutationResults_HbO,['clusterPerm_10_16_N_TLvsEL_HbO' '_nperm1K']);
writetable(ClusterPermutationResults_HbR,['clusterPerm_10_16_N_TLvsEL_HbR' '_nperm1K']); 


