% SMB_GLACIOCLIM.m: script to pre-process the SMB GLACIOCLIM measurements
% Author: Marin Kneib (marin.kneib@univ-grenoble-alpes.fr

clear
close all

% path to useful codes
addpath('../UsefulCodes/ImGRAFT')
addpath('../UsefulCodes/')

cd C:\Users\kneibm\Documents\CAIRN\Accu_Argentiere\argentiere_pleiades_smb\code

% folders
datafolder='../data/velocity/PLEIADES_ALPS/';
resultsfolder = '../output/smb_glacioclim/';

%% Calculate average SMB over 2012-2021 for each stake number

% options to load csvs
opts = delimitedTextImportOptions("NumVariables", 13);
opts.DataLines = [2, Inf];
opts.Delimiter = ";";
opts.VariableNames = ["profile_name", "Var2", "stake_number", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "x_lambert2e", "y_lambert2e", "altitude", "annual_smb"];
opts.SelectedVariableNames = ["profile_name", "stake_number", "x_lambert2e", "y_lambert2e", "altitude", "annual_smb"];
opts.VariableTypes = ["double", "string", "double", "string", "string", "string", "string", "string", "string", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["Var2", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var2", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "profile_name", "TrimNonNumeric", true);
opts = setvaropts(opts, "profile_name", "ThousandsSeparator", ",");

% Average measurements per stake location
SMB_abl = zeros(22,7);
SMB_acc = zeros(8,7);
SMB_tna = zeros(8,7);

All_abl = nan(22,10);
All_acc = nan(8,10);
All_tna = nan(8,10);

N_abl = zeros(22,1);
N_acc = zeros(8,1);
N_tna = zeros(8,1);

stakes_p4 = [1,2,3,4,5];
stakes_p5 = [1,2,3,4,5];
stakes_p7 = [1,2,3,4,5,6,7,8,9,10,11,12];
stakes_acc = [1,2,3,5,6,7,9,10];
stakes_tna = [1,2,3,4,5,6,7,8];

for yr = 2012:2021
    yr_str = string(yr);
    SMB_abl_yr = readtable(['..\data\smb_glacioclim\ablation_zone\Argentiere_Profils_2_4_5_7_annual_smb_abl_',yr_str{:},'.csv'], opts);
    SMB_acc_yr = readtable(['..\data\smb_glacioclim\accu_zone\Argentiere_annual_smb_accu_',yr_str{:},'.csv'], opts);
    SMB_tna_yr = readtable(['..\data\smb_glacioclim\amethystes_tournoir\Argentiere_Tour_Noir_annual_smb_abl_',yr_str{:},'.csv'], opts);
    idx = 0;
    for ss = stakes_p4
        pp = 4;
        idx = idx+1;
        meas = SMB_abl_yr.profile_name==pp & SMB_abl_yr.stake_number==ss;
        if sum(meas(:))>0
            N_abl(idx,1) = N_abl(idx,1)+1;
            SMB_abl(idx,1) = pp;
            SMB_abl(idx,2) = ss;
            SMB_abl(idx,3) = SMB_abl(idx,3)+mean(SMB_abl_yr.x_lambert2e(meas));
            SMB_abl(idx,4) = SMB_abl(idx,4)+mean(SMB_abl_yr.y_lambert2e(meas));
            SMB_abl(idx,5) = SMB_abl(idx,5)+mean(SMB_abl_yr.altitude(meas));
            SMB_abl(idx,6) = SMB_abl(idx,6)+mean(SMB_abl_yr.annual_smb(meas));
            All_abl(idx,yr-2011) = mean(SMB_abl_yr.annual_smb(meas));
        end
    end
    for ss = stakes_p5
        pp = 5;
        idx = idx+1;
        meas = SMB_abl_yr.profile_name==pp & SMB_abl_yr.stake_number==ss;
        if sum(meas(:))>0
            N_abl(idx,1) = N_abl(idx,1)+1;
            SMB_abl(idx,1) = pp;
            SMB_abl(idx,2) = ss;
            SMB_abl(idx,3) = SMB_abl(idx,3)+mean(SMB_abl_yr.x_lambert2e(meas));
            SMB_abl(idx,4) = SMB_abl(idx,4)+mean(SMB_abl_yr.y_lambert2e(meas));
            SMB_abl(idx,5) = SMB_abl(idx,5)+mean(SMB_abl_yr.altitude(meas));
            SMB_abl(idx,6) = SMB_abl(idx,6)+mean(SMB_abl_yr.annual_smb(meas));
            All_abl(idx,yr-2011) = mean(SMB_abl_yr.annual_smb(meas));
        end
    end
    for ss = stakes_p7
        pp = 7;
        idx = idx+1;
        meas = SMB_abl_yr.profile_name==pp & SMB_abl_yr.stake_number==ss;
        if sum(meas(:))>0
            N_abl(idx,1) = N_abl(idx,1)+1;
            SMB_abl(idx,1) = pp;
            SMB_abl(idx,2) = ss;
            SMB_abl(idx,3) = SMB_abl(idx,3)+mean(SMB_abl_yr.x_lambert2e(meas));
            SMB_abl(idx,4) = SMB_abl(idx,4)+mean(SMB_abl_yr.y_lambert2e(meas));
            SMB_abl(idx,5) = SMB_abl(idx,5)+mean(SMB_abl_yr.altitude(meas));
            SMB_abl(idx,6) = SMB_abl(idx,6)+mean(SMB_abl_yr.annual_smb(meas));
            All_abl(idx,yr-2011) = mean(SMB_abl_yr.annual_smb(meas));
        end
    end
    idx = 0;
    for ss = stakes_acc
        pp = NaN;
        idx = idx+1;
        meas = SMB_acc_yr.stake_number==ss;
        if sum(meas(:))>0
            N_acc(idx,1) = N_acc(idx,1)+1;
            SMB_acc(idx,1) = pp;
            SMB_acc(idx,2) = ss;
            SMB_acc(idx,3) = SMB_acc(idx,3)+mean(SMB_acc_yr.x_lambert2e(meas));
            SMB_acc(idx,4) = SMB_acc(idx,4)+mean(SMB_acc_yr.y_lambert2e(meas));
            SMB_acc(idx,5) = SMB_acc(idx,5)+mean(SMB_acc_yr.altitude(meas));
            SMB_acc(idx,6) = SMB_acc(idx,6)+mean(SMB_acc_yr.annual_smb(meas));
            All_acc(idx,yr-2011) = mean(SMB_acc_yr.annual_smb(meas));
        end
    end
    idx = 0;
    for ss = stakes_tna
        pp = NaN;
        idx = idx+1;
        meas = SMB_tna_yr.stake_number==ss;
        if sum(meas(:))>0
            N_tna(idx,1) = N_tna(idx,1)+1;
            SMB_tna(idx,1) = pp;
            SMB_tna(idx,2) = ss;
            SMB_tna(idx,3) = SMB_tna(idx,3)+mean(SMB_tna_yr.x_lambert2e(meas));
            SMB_tna(idx,4) = SMB_tna(idx,4)+mean(SMB_tna_yr.y_lambert2e(meas));
            SMB_tna(idx,5) = SMB_tna(idx,5)+mean(SMB_tna_yr.altitude(meas));
            SMB_tna(idx,6) = SMB_tna(idx,6)+mean(SMB_tna_yr.annual_smb(meas));
            All_tna(idx,yr-2011) = mean(SMB_tna_yr.annual_smb(meas));
        end
    end
end

% remove values where less than 50% of measurements and calculate temporal
% average
for ii = 1:size(N_abl,1)
    if N_abl(ii,1)<6
        SMB_abl(ii,:) = NaN;
    else
        for jj = 3:6
            SMB_abl(ii,jj) = SMB_abl(ii,jj)/N_abl(ii,1);
        end
        SMB_abl(ii,7) = nanstd(All_abl(ii,:));
    end
end
for ii = 1:size(N_acc,1)
    if N_acc(ii,1)<6
        SMB_acc(ii,:) = NaN;
    else
        for jj = 3:6
            SMB_acc(ii,jj) = SMB_acc(ii,jj)/N_acc(ii,1);
        end
        SMB_acc(ii,7) = nanstd(All_acc(ii,:));
    end
end
for ii = 1:size(N_tna,1)
    if N_tna(ii,1)<6
        SMB_tna(ii,:) = NaN;
    else
        for jj = 3:6
            SMB_tna(ii,jj) = SMB_tna(ii,jj)/N_tna(ii,1);
        end
        SMB_tna(ii,7) = nanstd(All_tna(ii,:));
    end
end

%% export temporal average of stake measurements as csv

writematrix(SMB_abl,[resultsfolder,'SMB_ABL_2012-2021_epsg27562.csv']);
writematrix(SMB_acc,[resultsfolder,'SMB_ACC_2012-2021_epsg27562.csv']);
writematrix(SMB_tna,[resultsfolder,'SMB_TNA_2012-2021_epsg27562.csv']);









