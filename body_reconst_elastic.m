% Reconstruct snake body using pre-processed data produced by
% Pre_processing.m

close all; clearvars;clc

%% Input

% The form recording marker order, orientation correction info etc.
infoForm = 'D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Data\Other\PreProcessingInfo.xlsx';
conditionLine = 2; % Which line in the excel form

% Preprocessed tag data file path, note that the unit of length data is mm
ppfilepath = 'D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Data\Pre-processed';
% Save path
savepath = 'D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Data\Reconstructed';

%% Reconstruction parameters
extra_dof = 3; % Extra DOF used for B-spline reconstruction
%tag correction, smoothing, and to body
th = 1.8; %height above body
s_m = 20.2; %mass of snake in grams
s_m_kg = .001*s_m; %convert to kg
rad = 4.5;%5.54;%radius of snake in mm
rad_m = .001*rad; %radius of snake in m
tot_lgth = 393.7; %length of snake in mm
tot_lgth_m = .001*tot_lgth; % length of snake in m
s_A = pi*(rad_m)^2; %cross sectional area of snake in m^2
s_V = tot_lgth_m*s_A; %total volume of snake m^3
den = s_m/s_V; %density of snake body in kg/m^3
grav = 9.8; %gravitational acceleration in m/s^2

Nseg = 701; % Number of small elements in each segment

%% Load information
% Marker number sequence (from head to tail)
[~,order,~] = xlsread(infoForm,'Sheet1',['E' num2str(conditionLine)]); % Tag number list
    order = str2num(order{1});
% Length adjustment for each segment
% If reconstruction doesn't look correct, manually adjust this value for
% better results
[~,lgth_adj,~] = xlsread(infoForm,'Sheet1',['F' num2str(conditionLine)]); % length adjustment proportion, negative if shorter
    lgth_adj = str2num(lgth_adj{1});

frameRange = []; % [starting ending] frame; empty if all frames need to be included

ppFiles = dir([ppfilepath,'\Trial-*.mat']);

%% Processing
for indFile = 1:length(ppFiles) % Each file
    tag_data_fullpath = [ppfilepath,'\',ppFiles(indFile).name]; % Path of this file
    
    % Check data
    ppData = load([ppfilepath,'\',ppFiles(indFile).name]);
    if length(order)~=size(ppData.XPlot,1)
        error('Tag number mismatch');
    end
    
    % Save as
    split = strsplit(ppFiles(indFile).name, '_pp');
    DATAID = split{1};
    save_name = [savepath,'\',DATAID, '_backbone.mat'];

    % RML = [radius (mm), mass (g), total length (mm)]
    RML = [rad, s_m, tot_lgth];

    reconstruct_body_elastic(tag_data_fullpath, save_name, RML, frameRange, lgth_adj);
end