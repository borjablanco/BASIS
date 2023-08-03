% -----------------------------------------
% Script for group results BASIS fNIRS ----
% Extract data from time window of interst-
% ----- Boxplots and correlation plots ----
% Borja Blanco - bb579@cam.ac.uk
% -----------------------------------------

% Initialize environment
clear; close all; clc

% Add required toolboxes/functions to path
addpath(genpath('/Users/borjablanco/GitHub/BASIS'))

% Load data from all participants (for each condition)
% Stored in variables with names as below

% group_N_hbo = [];
% group_N_hbr = [];
% group_V_hbo = [];
% group_V_hbr = [];
% group_S_hbo = [];
% group_S_hbr = [];

% Group data will be in format (Time x Channels x Participants)
% Time = 241
% Channels = 26
% Participants = 81
% At this point, channels, trials and subjects have already been excluded

% Number of participants in each group
ncont = 23;
nhfl = 58;

% Extract data from channels belonging to a specific cluster
cluster = [12,13,14,16,18];
tw = 141:201; % 10-16 seconds

% Extract data from control group
n_cont_hbo = mean(squeeze(nanmean(group_N_hbo(tw,cluster,1:ncont),2)),1);
n_cont_hbr = mean(squeeze(nanmean(group_N_hbr(tw,cluster,1:ncont),2)),1);
v_cont_hbo = mean(squeeze(nanmean(group_V_hbo(tw,cluster,1:ncont),2)),1);
v_cont_hbr = mean(squeeze(nanmean(group_V_hbr(tw,cluster,1:ncont),2)),1);
s_cont_hbo = mean(squeeze(nanmean(group_S_hbo(tw,cluster,1:ncont),2)),1);
s_cont_hbr = mean(squeeze(nanmean(group_S_hbr(tw,cluster,1:ncont),2)),1);

n_hfl_hbo = mean(squeeze(nanmean(group_N_hbo(tw,cluster,ncont+1:end),2)),1);
n_hfl_hbr = mean(squeeze(nanmean(group_N_hbr(tw,cluster,ncont+1:end),2)),1);
v_hfl_hbo = mean(squeeze(nanmean(group_V_hbo(tw,cluster,ncont+1:end),2)),1);
v_hfl_hbr = mean(squeeze(nanmean(group_V_hbr(tw,cluster,ncont+1:end),2)),1);
s_hfl_hbo = mean(squeeze(nanmean(group_S_hbo(tw,cluster,ncont+1:end),2)),1);
s_hfl_hbr = mean(squeeze(nanmean(group_S_hbr(tw,cluster,ncont+1:end),2)),1);

% ADHD only 
%adhd_idx = [];
adhd_sub = sub(adhd_idx);
n_adhd_hbo = mean(squeeze(nanmean(group_N_hbo(tw,cluster,adhd_idx),2)),1);
n_adhd_hbr = mean(squeeze(nanmean(group_N_hbr(tw,cluster,adhd_idx),2)),1);
v_adhd_hbo = mean(squeeze(nanmean(group_V_hbo(tw,cluster,adhd_idx),2)),1);
v_adhd_hbr = mean(squeeze(nanmean(group_V_hbr(tw,cluster,adhd_idx),2)),1);
s_adhd_hbo = mean(squeeze(nanmean(group_S_hbo(tw,cluster,adhd_idx),2)),1);
s_adhd_hbr = mean(squeeze(nanmean(group_S_hbr(tw,cluster,adhd_idx),2)),1);

% ASD+ADHD
%asdadhd_idx = [];
asdadhd_sub = sub(asdadhd_idx);
n_asdadhd_hbo = mean(squeeze(nanmean(group_N_hbo(tw,cluster,asdadhd_idx),2)),1);
n_asdadhd_hbr = mean(squeeze(nanmean(group_N_hbr(tw,cluster,asdadhd_idx),2)),1);
v_asdadhd_hbo = mean(squeeze(nanmean(group_V_hbo(tw,cluster,asdadhd_idx),2)),1);
v_asdadhd_hbr = mean(squeeze(nanmean(group_V_hbr(tw,cluster,asdadhd_idx),2)),1);
s_asdadhd_hbo = mean(squeeze(nanmean(group_S_hbo(tw,cluster,asdadhd_idx),2)),1);
s_asdadhd_hbr = mean(squeeze(nanmean(group_S_hbr(tw,cluster,asdadhd_idx),2)),1);

% ASD only
%asd_idx = 1:length(sub);
asd_idx([1:ncont, adhd_idx, asdadhd_idx]) = [];
asd_sub = sub(asd_idx);
n_asd_hbo = mean(squeeze(nanmean(group_N_hbo(tw,cluster,asd_idx),2)),1);
n_asd_hbr = mean(squeeze(nanmean(group_N_hbr(tw,cluster,asd_idx),2)),1);
v_asd_hbo = mean(squeeze(nanmean(group_V_hbo(tw,cluster,asd_idx),2)),1);
v_asd_hbr = mean(squeeze(nanmean(group_V_hbr(tw,cluster,asd_idx),2)),1);
s_asd_hbo = mean(squeeze(nanmean(group_S_hbo(tw,cluster,asd_idx),2)),1);
s_asd_hbr = mean(squeeze(nanmean(group_S_hbr(tw,cluster,asd_idx),2)),1);

% DATA FOR CONTRAST
vn_cont_hbo = v_cont_hbo - n_cont_hbo;
vn_cont_hbr = v_cont_hbr - n_cont_hbr;

vn_hfl_hbo = v_hfl_hbo - n_hfl_hbo;
vn_hfl_hbr = v_hfl_hbr - n_hfl_hbr;

vn_asd_hbo = v_asd_hbo - n_asd_hbo;
vn_asd_hbr = v_asd_hbr - n_asd_hbr;

vn_adhd_hbo = v_adhd_hbo - n_adhd_hbo;
vn_adhd_hbr = v_adhd_hbr - n_adhd_hbr;

vn_asdadhd_hbo = v_asdadhd_hbo - n_asdadhd_hbo;
vn_asdadhd_hbr = v_asdadhd_hbr - n_asdadhd_hbr;

% Create idxs for each group
cont_num = ones(ncont,1);
asd_num = 2*ones(size(asd_sub,1),1);
adhd_num = 3*ones(size(adhd_sub,1),1);
asdadhd_num = 4*ones(size(asdadhd_sub,1),1);
group_num = [cont_num; asd_num; adhd_num; asdadhd_num];

% Store data for plots 
stats_group_s_hbo = [s_cont_hbo'; s_asd_hbo'; s_adhd_hbo'; s_asdadhd_hbo'];
stats_group_s_hbr = [s_cont_hbr'; s_asd_hbr'; s_adhd_hbr'; s_asdadhd_hbr'];


% 2 - Create boxplot figure for HbO and HbR
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);

subplot(5,2,[1,3,5])
bb_boxplots(group_num, stats_group_s_hbo, cont_num, asd_num, adhd_num, asdadhd_num, asd_code, asdadhd_code)
title('HbO'); 
ylim([-1 1.5])

subplot(5,2,[2,4,6])
bb_boxplots(group_num, stats_group_s_hbr, cont_num, asd_num, adhd_num, asdadhd_num, asd_code, asdadhd_code)
title('HbR');
ylim([-0.75 0.5])

% 3 - Correlation plots

% Load table with data
%basis = readtable('');

% Extract info
id = basis.SIBSID;
group = basis.Group;
cluster = basis.clusterVN_HbO;
behav = basis.SRS2_Total;
asddiag = basis.ASDDiag;

% Group variables
data = [behav, cluster, group, asddiag];

% Sort by group (for plots)
[~, idx] = sort(data(:,3),'ascend');
sortdata = data(idx,:);

% Remove NaN
plotdata = sortdata;
nan_behav = isnan(plotdata(:,1));
nan_cluster = isnan(plotdata(:,2));
nan_idx = or(nan_behav, nan_cluster); 
plotdata(find(nan_idx),:) = [];

% Assign idx to groups
[~,b] = unique(plotdata(:,3));
tl = 1:b(2)-1;
asd = b(2):b(3)-1;
adhd = b(3):b(4)-1;
asdadhd = b(4):size(plotdata,1);

A = plotdata(:,1);
B = plotdata(:,2);

fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
subplot(2,2,1)
Po = polyfit(A,B,1);     % 1st order polynomial (linear)
Yfit = Po(1)*A+Po(2);    % y-values of the fitted line
plot(A,Yfit,'r', 'Linewidth',2);       % Fitted line (color: red, type: solid)
hold on
scatter(A(plotdata(:,4)==1,:),B(plotdata(:,4)==1), 300, 'Marker', '+', 'MarkerEdgeColor', 'r', 'Linewidth', 2);
scatter(A(tl),B(tl),80, 'filled', 'black', 'MarkerEdgeColor', 'black', 'Linewidth',1)
scatter(A(asd),B(asd),80,'filled', 'b', 'Marker', 'square', 'MarkerEdgeColor', 'b', 'Linewidth',1)
scatter(A(adhd),B(adhd),80,'MarkerFaceColor', [1, 0.7, 0], 'Marker', 'diamond', 'MarkerEdgeColor', [1, 0.7, 0], 'Linewidth',1)
scatter(A(asdadhd),B(asdadhd),80,'MarkerFaceColor', 'g', 'Marker', '*',  'MarkerEdgeColor', 'g', 'Linewidth',1)
set(gca, 'Fontsize', 24)
title (['SRS2 Total - Cluster VN HbO, r = ' num2str(corr(A,B, 'Type', 'spearman'))], 'FontSize',20)
[r,p] = corr(A,B, 'Type', 'spearman');
box off
ylim([-1.6 1.6])

% Repeat for each behav measure / cluster

