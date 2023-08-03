% -----------------------------------------
% Script for group results BASIS fNIRS ----
% Borja Blanco - bb579@cam.ac.uk
% -----------------------------------------
% Initialize environment
clear; close all; clc

% Add required toolboxes/functions to path
addpath(genpath('/Users/borjablanco/GitHub/BASIS'))
addpath(genpath('/Users/borjablanco/Documents/MATLAB/homer2'))

% Paths to preprocessed data
path_prep_control = '/Volumes/CAM_data/STAARS/BB_analysis/analysis_072021/control/preprocessed';
path_prep_highRisk = '/Volumes/CAM_data/STAARS/BB_analysis/analysis_072021/highRisk/preprocessed';

% Paths for storing figures
path_figures = '/Volumes/CAM_data/STAARS/BB_analysis/figures';

% List participants 
cd(path_prep_control)
sub_cont = dir('*prepro.mat');

cd(path_prep_highRisk)
sub_highRisk = dir('*prepro.mat');

% CONTROL GROUP
prep_path = path_prep_control;
cd(prep_path)

% EXCLUDE DATA FROM THOSE CHANNELS MARKED AS NOISY DURING PREPROCESSING - data.bad_links
nch = 26;
bad_channels_c = zeros(nch, length(sub_cont));
ntrials_c = zeros(length(sub_cont), 3);

% Initialize variables to store data
group_S_hbo_c = [];
group_S_hbr_c = [];
group_N_hbo_c = [];
group_N_hbr_c = [];
group_V_hbo_c = [];
group_V_hbr_c = [];

% For each participant
for nsub = 1:length(sub_cont)
    
    % Load data
    data = load(sub_cont(nsub).name);
    data = data.data;
        
    % Find idx of bad channels (only 1 wavelength)
    idx_ch = data.bad_links(1:length(data.bad_links)/2);
        
    % Exclude bad channels for group analysis
    hbo_s = data.hbo_S; hbo_s(:, idx_ch) = NaN;
    hbr_s = data.hbr_S; hbr_s(:, idx_ch) = NaN;
    hbo_n = data.hbo_N; hbo_n(:, idx_ch) = NaN;
    hbr_n = data.hbr_N; hbr_n(:, idx_ch) = NaN;
    hbo_v = data.hbo_V; hbo_v(:, idx_ch) = NaN;
    hbr_v = data.hbr_V; hbr_v(:, idx_ch) = NaN;
    
    bad_channels_c (idx_ch, nsub) = 1;
    
    % Load estimated HRF
    group_S_hbo_c = cat(3, group_S_hbo_c, hbo_s);
    group_S_hbr_c = cat(3, group_S_hbr_c, hbr_s);
    group_N_hbo_c = cat(3, group_N_hbo_c, hbo_n);
    group_N_hbr_c = cat(3, group_N_hbr_c, hbr_n);
    group_V_hbo_c = cat(3, group_V_hbo_c, hbo_v);
    group_V_hbr_c = cat(3, group_V_hbr_c, hbr_v);
    
    % Order as V, N, S
    ntrials_c(nsub,:) = data.nTrials([4,3,2]);
    
end

ch_disc_c = sum(bad_channels_c)';

% ELEVATED LIKELIHOOD GROUP
prep_path = path_prep_highRisk;
cd(prep_path)

% EXCLUDE DATA FROM THOSE CHANNELS MARKED AS NOISY DURING PREPROCESSING - data.bad_links
bad_channels_e = zeros(nch, length(sub_highRisk));
ntrials_e = zeros(length(sub_highRisk), 3);

% Initialize variables to store data
group_S_hbo_e = [];
group_S_hbr_e = [];
group_N_hbo_e = [];
group_N_hbr_e = [];
group_V_hbo_e = [];
group_V_hbr_e = [];

% For each participant
for nsub = 1:length(sub_highRisk)
    
    % Load data
    data = load (sub_highRisk(nsub).name);
    data = data.data;
        
    % Find idx of bad channels (only 1 wavelength)
    idx_ch = data.bad_links(1:length(data.bad_links)/2);
        
    % Exclude bad channels for group analysis
    hbo_s = data.hbo_S; hbo_s(:, idx_ch) = NaN;
    hbr_s = data.hbr_S; hbr_s(:, idx_ch) = NaN;
    hbo_n = data.hbo_N; hbo_n(:, idx_ch) = NaN;
    hbr_n = data.hbr_N; hbr_n(:, idx_ch) = NaN;
    hbo_v = data.hbo_V; hbo_v(:, idx_ch) = NaN;
    hbr_v = data.hbr_V; hbr_v(:, idx_ch) = NaN;
    
    bad_channels_e (idx_ch, nsub) = 1;
    
    % Load estimated HRF
    group_S_hbo_e = cat(3, group_S_hbo_e, hbo_s);
    group_S_hbr_e = cat(3, group_S_hbr_e, hbr_s);
    group_N_hbo_e = cat(3, group_N_hbo_e, hbo_n);
    group_N_hbr_e = cat(3, group_N_hbr_e, hbr_n);
    group_V_hbo_e = cat(3, group_V_hbo_e, hbo_v);
    group_V_hbr_e = cat(3, group_V_hbr_e, hbr_v);
       
    ntrials_e(nsub,:) = data.nTrials([4,3,2]);
    
end

ch_disc_e = sum(bad_channels_e)';

% Concatenate control and elevated likelihood groups
group_S_hbo = cat(3, group_S_hbo_c, group_S_hbo_e);
group_S_hbr = cat(3, group_S_hbr_c, group_S_hbr_e);
group_N_hbo = cat(3, group_N_hbo_c, group_N_hbo_e);
group_N_hbr = cat(3, group_N_hbr_c, group_N_hbr_e);
group_V_hbo = cat(3, group_V_hbo_c, group_V_hbo_e);
group_V_hbr = cat(3, group_V_hbr_c, group_V_hbr_e);


% ########## DISCARD PARTICIPANTS WITH A LOW NUMBER OF CHANNELS INCLUDED
ch_discard = find(sum(bad_channels)'>20);

% ########## DISCARD PARTICIPANTS WITH A LOW NUMBER OF TRIALS INCLUDED
tr_discard = find(ntrials(1,:)<3 | ntrials(2,:)<3 | ntrials(3,:)<3);

% exclude participants
exc_idx = unique([ch_discard, tr_discard]);

% Exclude participants based on n trials and n channels
group_S_hbo(:,:,exc_idx) = [];
group_S_hbr(:,:,exc_idx) = [];
group_N_hbo(:,:,exc_idx) = [];
group_N_hbr(:,:,exc_idx) = [];
group_V_hbo(:,:,exc_idx) = [];
group_V_hbr(:,:,exc_idx) = [];

% Compute group-averaged HRF response
avg_S_hbo = nanmean(group_S_hbo,3);
avg_S_hbr = nanmean(group_S_hbr,3);
avg_N_hbo = nanmean(group_N_hbo,3);
avg_N_hbr = nanmean(group_N_hbr,3);
avg_V_hbo = nanmean(group_V_hbo,3);
avg_V_hbr = nanmean(group_V_hbr,3);

% Compute standard deviation for plots
std_S_hbo = std(group_S_hbo,[],3, 'omitnan');
std_S_hbr = std(group_S_hbr,[],3, 'omitnan');
std_N_hbo = std(group_N_hbo,[],3, 'omitnan');
std_N_hbr = std(group_N_hbr,[],3, 'omitnan');
std_V_hbo = std(group_V_hbo,[],3, 'omitnan');
std_V_hbr = std(group_V_hbr,[],3, 'omitnan');

% PLOT AVERAGE RESPONSES AT THE GROUP LEVEL
trange = [-4 20];
sf = 10;
time = (trange(1)*sf:ceil(trange(2)*sf))/sf;

% Locations for plots based on SD file
locations = [2 8 10 16 24 22 30 36 38 44 52 6 14 12 26 ...
    20 34 28 42 40 54 48 62 56 50 58];

% VOCAL CONDITION
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
for i = 1:size(data.OD,2)/2
    sphand = subplot(9,7,locations(i));
    pos1 = get(sphand, 'Position'); % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 0.02];
    set(sphand, 'Position',new_pos1 ) % set new position of current sub - plot
        
    s = shadedErrorBar(time, avg_V_hbo (:,i), std_V_hbo(:,i), 'lineprops', '-r' ,'patchsaturation', 0.1);
    set(s.mainLine,'LineWidth',1.5)
    hold on
    s = shadedErrorBar(time, avg_V_hbr (:,i), std_V_hbr(:,i), 'lineprops', '-blue','patchsaturation', 0.1);
    set(s.mainLine,'LineWidth',1.5)
    
    area([0, 8], [2 2], 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'LineStyle', 'none')
    
    title (['ch ' num2str(i)], 'fontsize', 16)
    set(gca, 'fontsize', 16)
    ylim([-1.2 1.7]);
    xlim([-4 20])
    box off
    
    if i == 21
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
    if i == 25
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
    
end

% NON-VOCAL CONDITION
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
for i = 1:size(data.OD,2)/2
    sphand = subplot(9,7,locations(i));
    pos1 = get(sphand, 'Position'); % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 0.02];
    set(sphand, 'Position',new_pos1 ) % set new position of current sub - plot
        
    s = shadedErrorBar(time, avg_N_hbo (:,i), std_N_hbo(:,i), 'lineprops', '-m' ,'patchsaturation', 0.1);
    set(s.mainLine,'LineWidth',1.5)
    hold on
    s = shadedErrorBar(time, avg_N_hbr (:,i), std_N_hbr(:,i), 'lineprops', '-c','patchsaturation', 0.1);
    set(s.mainLine,'LineWidth',1.5)
    
    area([0, 8], [2 2], 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'LineStyle', 'none')
    
    title (['ch ' num2str(i)], 'fontsize', 16)
    set(gca, 'fontsize', 16)
    ylim([-1.2 1.7]);
    xlim([-4 20])
    box off
    
    if i == 21
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
    if i == 25
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
    
end

% VISUAL SOCIAL CONDITION
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
for i = 1:size(data.OD,2)/2
    sphand = subplot(9,7,locations(i));
    pos1 = get(sphand, 'Position'); % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 0.02];
    set(sphand, 'Position',new_pos1 ) % set new position of current sub - plot
        
    s = shadedErrorBar(time, avg_S_hbo (:,i), std_S_hbo(:,i), 'lineprops', '-g' ,'patchsaturation', 0.1);
    set(s.mainLine,'LineWidth',1.5)
    hold on
    s = shadedErrorBar(time, avg_S_hbr (:,i), std_S_hbr(:,i), 'lineprops', {'-w','markerfacecolor',[1, 0.7, 0]}, 'patchsaturation', 0.1);
    set(s.mainLine,'LineWidth',1.5, 'color', [1, 0.7, 0])
    set(s.edge,'LineWidth',0.5, 'color', [1, 0.7, 0])
    set(s.patch, 'FaceColor', [1, 0.7, 0]);

    
    area([0, 8], [2 2], 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'LineStyle', 'none')
    %area([0, 8], [-2 -2], 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'LineStyle', 'none')

    title (['ch ' num2str(i)], 'fontsize', 16)
    set(gca, 'fontsize', 16)
    ylim([-1 1.4]);
    xlim([-4 20])
    box off
    
    if i == 21
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
    if i == 25
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
    
end

% V vs NV CONDITIONS
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
for i = 1:size(data.OD,2)/2
    hold on
    sphand = subplot(9,7,locations(i));
    pos1 = get(sphand, 'Position'); % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 0.02];
    set(sphand, 'Position',new_pos1 ) % set new position of current sub - plot
            
    % FOR GAUSSIANS
    plot(time, avg_N_hbo(:,i), 'm', 'linewidth', 1.5, 'LineStyle', '-')
    hold on
    plot(time,avg_N_hbr(:,i),'c', 'linewidth', 1.5, 'LineStyle', '-')
    plot(time, avg_V_hbo(:,i), 'r', 'linewidth', 1.5, 'LineStyle', '-')
    plot(time, avg_V_hbr(:, i), 'blue', 'linewidth', 1.5, 'LineStyle', '-')
    
    area([0, 8], [2 2], 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'LineStyle', 'none')
    
    title (['ch ' num2str(i)], 'fontsize', 16)
    set(gca, 'fontsize', 16)
    ylim([-0.6 1.2]);
    xlim([-4 20])
    box off
    
    if i == 21
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
    if i == 25
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
end

% PLOTS FOR EACH GROUP
% CONTROL
% Compute group-averaged HRF response
avg_S_hbo_c = nanmean(group_S_hbo_c,3);
avg_S_hbr_c = nanmean(group_S_hbr_c,3);
avg_N_hbo_c = nanmean(group_N_hbo_c,3);
avg_N_hbr_c = nanmean(group_N_hbr_c,3);
avg_V_hbo_c = nanmean(group_V_hbo_c,3);
avg_V_hbr_c = nanmean(group_V_hbr_c,3);

% Compute standard deviation for plots
std_S_hbo_c = std(group_S_hbo_c,[],3, 'omitnan');
std_S_hbr_c = std(group_S_hbr_c,[],3, 'omitnan');
std_N_hbo_c = std(group_N_hbo_c,[],3, 'omitnan');
std_N_hbr_c = std(group_N_hbr_c,[],3, 'omitnan');
std_V_hbo_c = std(group_V_hbo_c,[],3, 'omitnan');
std_V_hbr_c = std(group_V_hbr_c,[],3, 'omitnan');

% ELEVATED LIKELIHOOD
% Compute group-averaged HRF response
avg_S_hbo_e = nanmean(group_S_hbo_e,3);
avg_S_hbr_e = nanmean(group_S_hbr_e,3);
avg_N_hbo_e = nanmean(group_N_hbo_e,3);
avg_N_hbr_e = nanmean(group_N_hbr_e,3);
avg_V_hbo_e = nanmean(group_V_hbo_e,3);
avg_V_hbr_e = nanmean(group_V_hbr_e,3);

% Compute standard deviation for plots
std_S_hbo_e = std(group_S_hbo_e,[],3, 'omitnan');
std_S_hbr_e = std(group_S_hbr_e,[],3, 'omitnan');
std_N_hbo_e = std(group_N_hbo_e,[],3, 'omitnan');
std_N_hbr_e = std(group_N_hbr_e,[],3, 'omitnan');
std_V_hbo_e = std(group_V_hbo_e,[],3, 'omitnan');
std_V_hbr_e = std(group_V_hbr_e,[],3, 'omitnan');

% V vs NV CONDITIONS
% And S CONDITION
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
for i = 1:size(data.OD,2)/2
    hold on
    sphand = subplot(9,7,locations(i));
    pos1 = get(sphand, 'Position'); % gives the position of current sub-plot
    new_pos1 = pos1 +[0 0 0 0.02];
    set(sphand, 'Position',new_pos1 ) % set new position of current sub - plot
            
    % Plot group averaged hemodynamic responses
    plot(time, avg_N_hbo_c(:,i), 'm', 'linewidth', 1.5, 'LineStyle', '-')
    hold on
    plot(time, avg_N_hbr_c(:,i),'c', 'linewidth', 1.5, 'LineStyle', '-')
    plot(time, avg_V_hbo_c(:,i), 'r', 'linewidth', 1.5, 'LineStyle', '-')
    plot(time, avg_V_hbr_c(:, i), 'blue', 'linewidth', 1.5, 'LineStyle', '-')
    
    plot(time, avg_N_hbo_e(:,i), 'm', 'linewidth', 1.5, 'LineStyle', '-.')
    hold on
    plot(time, avg_N_hbr_e(:,i),'c', 'linewidth', 1.5, 'LineStyle', '-.')
    plot(time, avg_V_hbo_e(:,i), 'r', 'linewidth', 1.5, 'LineStyle', '-.')
    plot(time, avg_V_hbr_e(:, i), 'blue', 'linewidth', 1.5, 'LineStyle', '-.')

    % FOR SILENCE CONDITION
    %plot(time, avg_S_hbo_c(:,i), 'g', 'linewidth', 1.5, 'LineStyle', '-')
    %hold on
    %plot(time,avg_S_hbr_c(:,i),'color', [1, 0.7, 0], 'linewidth', 1.5, 'LineStyle', '-')
    %plot(time, avg_S_hbo_e(:,i), 'g', 'linewidth', 1.5, 'LineStyle', '-.')
    %plot(time,avg_S_hbr_e(:,i),'color', [1, 0.7, 0], 'linewidth', 1.5, 'LineStyle', '-.')

    area([0, 8], [2 2], 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'LineStyle', 'none')
    
    title (['ch ' num2str(i)], 'fontsize', 16)
    set(gca, 'fontsize', 16)
    ylim([-0.5 1.1]);
    xlim([-4 20])
    box off
    
    if i == 21
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
    if i == 25
        xlabel('Time (seconds)')
        ylabel('Conc. change (\muM)')
    end
end


