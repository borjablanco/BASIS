function[obs_val_HbO, perm_vals_HbO, obs_val_HbR, perm_vals_HbR] = ClusterPermutationAnalysis_Basis_Comp(HbO, HbR, ROI_channels, NPerms, nTL, nEL)

% References:
% Maris & Oostenveld, Journal of Neuroscience Methods 2007
% Ferry et al Developmental Science 2016
% Benavides-Varela & Gervain DCN 2017
% Abboub et al Brain & Lanfuage 2016
% Implemented from the algorithm suggested by Cohen, M. Chapter 15, "Matlab for Cognitive Neuroscientists 2017
% 
% Anna Blasi Sept 2018HbR
% Borja Blanco Jan 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

% Cluster information:
% ROI_channels = [14,15,16,17,18];
% obs_val_HbO2 = 15.3217;
%ROI_channels = [14,15,16];
%obs_val_HbO2 = 10.4997;


[Nbabies, Nchannels] = size(HbO);

obs_val_HbO = 0;
obs_val_HbR = 0;

for tr = ROI_channels(1:end)
    [~,~,~,STATS_HbO] = ttest2(HbO(1:nTL, tr), HbO(nTL+1:end, tr));% Only channels with valid data
    obs_val_HbO = obs_val_HbO + STATS_HbO.tstat; %t-values;
    
    [~,~,~,STATS_HbR] = ttest2(HbR(1:nTL, tr), HbR(nTL+1:end, tr));% Only channels with valid data    
    obs_val_HbR = obs_val_HbR + STATS_HbR.tstat; %t-values;
end

% Separate HbO2 and HbR data of each participant 
%HbO = TW(:,2:2:end);
%HbR = TW(:,3:2:end);

% Generate Nperms randomized shuffles of the data set:
%NPerms = 10000;

% Use subsample of babies in EL (n=61) to make it equal to TL (n=25)
%HbO_TL = HbO(1:nTL,:);
%HbR_TL = HbR(1:nTL,:);
%HbO_EL = HbO(nTL+1:end,:);
%HbR_EL = HbR(nTL+1:end,:);

%sub_idx = randperm(nEL, nTL);
%HbO_EL = HbO_EL(sub_idx,:);
%HbR_EL = HbR_EL(sub_idx,:);

%HbO = [HbO_TL; HbO_EL];
%HbR = [HbR_TL; HbR_EL];
%Nbabies_sub = size(HbO,1);
% -----------------
perm_vals_HbO(1:NPerms) = zeros;
perm_vals_HbR(1:NPerms) = zeros;

for permi = 1:NPerms
        
    % Shuffle data
    fbabies = randperm(Nbabies); % Permute participant positions
    fchannels = randperm(Nchannels); % for all channels
  
    HbO_Tstat_channel(1:Nchannels) = zeros;
    HbR_Tstat_channel(1:Nchannels) = zeros;

    for k = 1:Nchannels

        ch = fchannels(k);

        % TTest stats
        [~,~,~,STATS_HbO] = ttest2(HbO(fbabies(1:nTL)', ch), HbO(fbabies(nTL+1:end)', ch));% Only channels with valid data
        HbO_Tstat_channel(k) = STATS_HbO.tstat;
        [~,~,~,STATS_HbR] = ttest2(HbR(fbabies(1:nTL)', ch), HbR(fbabies(nTL+1:end)', ch));% Only channels with valid data
        HbR_Tstat_channel(k) = STATS_HbR.tstat;

    end
    
    % Output T
    perm_vals_HbO(permi) = sum(HbO_Tstat_channel(ROI_channels)); % sum of all t-values from channels in the cluster
    perm_vals_HbR(permi) = sum(HbR_Tstat_channel(ROI_channels)); % sum of all t-values from channels in the cluster
     
     
     clear STATS_HbO TF_HbO STATS_HbR TF_HbR
     
end


% Plot histograms of the calculated distribution (one per HbO2 and one per
% HbO

%[N,X] = hist(perm_vals_HbO2,100);
% figure; hist(perm_vals_HbO2,100);
% hold on;
% plot([obs_val_HbO2 obs_val_HbO2], get(gca,'ylim')); hold off;
%P_HbO2 = sum(N(find(X>obs_val_HbO2))) / sum(N);






% 
% hist(perm_vals_HbR,100)
% hold on;
% plot([obs_valHbR obs_valHbR], get(gca,'ylim')); hold off;
