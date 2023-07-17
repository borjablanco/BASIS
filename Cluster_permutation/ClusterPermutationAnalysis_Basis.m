function[obs_val_HbO, P_HbO, obs_val_HbR, P_HbR] = ClusterPermutationAnalysis_Basis(HbO, HbR, ROI_channels, NPerms)

% References:
% Maris & Oostenveld, Journal of Neuroscience Methods 2007
% Ferry et al Developmental Science 2016
% Benavides-Varela & Gervain DCN 2017
% Abboub et al Brain & Lanfuage 2016
% Implemented from the algorithm suggested by Cohen, M. Chapter 15, "Matlab for Cognitive Neuroscientists 2017
% 
% Anna Blasi Sept 2018
% Borja Blanco Jan 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

% Cluster information:
% ROI_channels = [14,15,16,17,18];
% obs_val_HbO2 = 15.3217;
%ROI_channels = [14,15,16];
%obs_val_HbO2 = 10.4997;

% t-test parameters:
M = 0;         % mean value
ALPHA = 0.05;  % significance level
TAIL = 0;      % mean is different from M 

[Nbabies, Nchannels] = size(HbO);

obs_val_HbO = 0;
obs_val_HbR = 0;

for tr = ROI_channels(1:end)
    [~,~,~,STATS_HbO] = ttest(HbO(:, tr), M, ALPHA, TAIL);% Only channels with valid data
    obs_val_HbO = obs_val_HbO + STATS_HbO.tstat; %t-values;
    
    [~,~,~,STATS_HbR] = ttest(HbR(:, tr), M, ALPHA, TAIL);% Only channels with valid data    
    obs_val_HbR = obs_val_HbR + STATS_HbR.tstat; %t-values;
end

% Initialize variables to store perm data
perm_vals_HbO(1:NPerms) = zeros;
perm_vals_HbR(1:NPerms) = zeros;

% For each permutation
for permi = 1:NPerms
        
    % Shuffle data
    fchannels = randperm(Nchannels); % for all channels
  
    HbO_Tstat_channel(1:Nchannels) = zeros;
    HbR_Tstat_channel(1:Nchannels) = zeros;

    for k = 1:Nchannels

        ch = fchannels(k);

        % TTest stats
        [~,~,~,STATS_HbO] = ttest(HbO(:, ch), M, ALPHA, TAIL);% Only channels with valid data
        HbO_Tstat_channel(k) = STATS_HbO.tstat;
        [~,~,~,STATS_HbR] = ttest(HbR(:, ch), M, ALPHA, TAIL);% Only channels with valid data
        HbR_Tstat_channel(k) = STATS_HbR.tstat;

    end
    
    % Output T
    perm_vals_HbO(permi) = sum(HbO_Tstat_channel(ROI_channels)); % sum of all t-values from channels in the cluster
    perm_vals_HbR(permi) = sum(HbR_Tstat_channel(ROI_channels)); % sum of all t-values from channels in the cluster
          
    clear STATS_HbO STATS_HbR
     
end


% Plot histograms of the calculated distribution (one per HbO and one per
% HbO

%[N,X] = histogram(abs(perm_vals_HbO),100);
%figure; histogram(abs(perm_vals_HbO),100);
%hold on;
%plot([abs(obs_val_HbO) abs(obs_val_HbO)], get(gca,'ylim')); hold off;
%P_HbO = sum(N(find(X>obs_val_HbO))) / sum(N);


%BB - p-val:
Tdist_HbO = sort(abs(perm_vals_HbO), 'ascend');
%P_HbO = (length(find(Tdist_HbO>=obs_val_HbO)))/NPerms;
P_HbO = (length(find(Tdist_HbO>=abs(obs_val_HbO))))/NPerms; % for 2 sided tests

Tdist_HbR = sort(abs(perm_vals_HbR), 'ascend');
%P_HbR = (length(find(Tdist_HbR<=obs_val_HbR)))/NPerms;
P_HbR = (length(find(Tdist_HbR>=abs(obs_val_HbR))))/NPerms; % for 2 sided tests


% hist(perm_vals_HbR,100)
% hold on;
% plot([obs_valHbR obs_valHbR], get(gca,'ylim')); hold off;
