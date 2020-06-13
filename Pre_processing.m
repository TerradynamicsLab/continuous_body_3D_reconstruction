%% Pre-process tracked marker position and orientation data for input into reconstruction algorithm

% In the tracked marker data file, there are 6 variables: 
%   XPlot, YPlot, ZPlot, RollPlot, PitchPlot, YawPlot
% Each variable is a m x n matrix, where m is the number of markers, n is
%   the number of video frames
% XPlot, YPlot, and ZPlot record x/y/z positions (unit: mm)
% RollPlot, PitchPlot, and YawPlot record ZYX Euler angles (unit: radian)

% If markers were tracked and recorded in other formats, they need to be
% converted to this format, or the corresponding sections in the codes need
% to be modified.

% Housekeeping
close all; clearvars; clc;

%% Input parameters

% Location of tapering data file specific to individual
taper_dat = ['D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\',...
    'To release with sample data\Data\Other\tapering.mat'];

% Folder where you want to save your data
save_folder = 'D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Data\Pre-processed';

% The form recording marker order, orientation correction info etc.
infoForm = 'D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Data\Other\PreProcessingInfo.xlsx';
conditionLine = 2; % Which line in the excel form
         
% Manual yaw correction value, nan if no manual correction
yaw_c = nan;
% Markers to ignore during reconstruction
kill_tags = []; 
% Starting and ending frame
st_ed = [];

th = 2.54; % Height of marker measured from top of the snake (mm)
rad = 4.5; % Radius of snake (mm)

%% Begin pre-processing

%% Load data
% Folder containing tracking data files from different trials to pre-process
[~,data_folder,~] = xlsread(infoForm,'Sheet1',['B' num2str(conditionLine)]); 
    data_folder = data_folder{1};    
% Data file of the correction reference trial
[~,rep_trial_dat,~] = xlsread(infoForm,'Sheet1',['C' num2str(conditionLine)]); 
    rep_trial_dat = rep_trial_dat{1};            
% Frames used for reference
[~,frc,~] = xlsread(infoForm,'Sheet1',['D' num2str(conditionLine)]); 
    frc = str2num(frc{1});
% Marker number sequence (from head to tail)
[~,order,~] = xlsread(infoForm,'Sheet1',['E' num2str(conditionLine)]); 
    order = str2num(order{1});

if length(yaw_c)==1 && isnan(yaw_c)
    yaw_c = nan(length(order),1);
end

% Load tapering data
load(taper_dat);

% Load correction reference trial
load(rep_trial_dat); % XPlot, Y.., Z.., Roll.., Pitch.., Yaw..

%% Remove tags listed in kill_tags

% Generate a string in the file name to represent killed tags
if isempty(kill_tags)
    kill_str = '';
else
    kill_str = '';
    for i = 1:length(kill_tags)
        kill_str = [kill_str, num2str(kill_tags(i)), '_'];
    end   
end
num_tags = length(setdiff(order, kill_tags));

% Indices of tags kept
count = 0;
for i = kill_tags
    count = count + 1; 
    kill_ind(count) = find(order == i); % Index of tags to remove (number from head to tail)
end
if ~isempty(kill_tags)
    new_tag_ind = setdiff(1:length(order), kill_ind);
else
    new_tag_ind = 1:length(order);
end

%% Calculate tag orientation correction and spacing

% Tag orientation correction
count = 0;
for j = new_tag_ind
    count = count + 1;
    
    if ~isnan(yaw_c(j)) % With manual yaw correction
        yR = axis2rot([0, 0, yaw_c(j)]); % Yaw correction term
        if isnan(frc(j)) % With manual yaw correction; No reference frame
            Rcrt{count} = yR; % Correct based on manual estimation
        else % With manual estimation input, and a reference frame
            % Apply yaw correction, then correct roll and pitch
            RIfrc = EULERZYX_int([RollPlot(j, frc(j)), PitchPlot(j, frc(j)), YawPlot(j, frc(j))])*yR;
            anglels = EULERZYXINV_int(RIfrc);
            RIfrcyaw = EULERZYX_int([0, 0, anglels(3)]);
            Rcrt{count} = transpose(RIfrc)*RIfrcyaw;
        end
    elseif ~isnan(frc(j)) % No manual yaw correction; With reference frame
        % Correct roll and pitch
         RIfrc = EULERZYX_int([RollPlot(j, frc(j)), PitchPlot(j, frc(j)), YawPlot(j, frc(j))]);
         RIfrcyaw = EULERZYX_int([0, 0, YawPlot(j, frc(j))]);
         Rcrt{count} = transpose(RIfrc)*RIfrcyaw;      
    else % No manual yaw correction; No reference frame
        Rcrt{count} = eye(3); % No correction
    end

end       
clear XPlot YPlot ZPlot RollPlot PitchPlot YawPlot

% Compute segment length using max tag distances
cd(data_folder)
files = dir(['*.mat']);
all_lengths = zeros(length(files), length(order) - length(kill_tags) - 1);
al_sub_min = zeros(length(files)-1, length(order) - length(kill_tags) - 1);
for k = 1:length(files) % Each trial
    load(files(k).name);    
    % Remove tags
    count = 0;
    for j = new_tag_ind % Each tag kept
        count = count + 1;      
        XPlot_new(count, :) = XPlot(j, :);
        YPlot_new(count, :) = YPlot(j, :);
        ZPlot_new(count, :) = ZPlot(j, :);
        RollPlot_new(count, :) = RollPlot(j, :);
        PitchPlot_new(count, :) = PitchPlot(j, :);
        YawPlot_new(count, :) = YawPlot(j, :);
    end
    XPlot = XPlot_new;
    YPlot = YPlot_new;
    ZPlot = ZPlot_new;
    RollPlot = RollPlot_new;
    PitchPlot = PitchPlot_new;
    YawPlot = YawPlot_new;
    % Find tag spacing by finding maximum distance in all frames
    for indSec = 1:size(XPlot, 1) - 1 % Each section between tags
        tag1 = indSec;
        tag2 = indSec + 1;
        for ind = 1:length(XPlot)
            if norm([XPlot(tag1, ind) - XPlot(tag2, ind),...
                     YPlot(tag1, ind) - YPlot(tag2, ind),...
                     ZPlot(tag1, ind) - ZPlot(tag2, ind)]) ...
                    > all_lengths(k, indSec)           
                all_lengths(k, indSec) = norm([XPlot(tag1, ind) - XPlot(tag2, ind), YPlot(tag1, ind) - YPlot(tag2, ind), ZPlot(tag1, ind) - ZPlot(tag2, ind)]);
            end
        end
    end    
    clear XPlot YPlot ZPlot RollPlot PitchPlot YawPlot XPlot_new YPlot_new ZPlot_new RollPlot_new PitchPlot_new YawPlot_new
end
% Remove the min values among all trials for each section
for i = 1:size(all_lengths, 2)
    al_sub_min(:, i) = setdiff(all_lengths(:, i), [min(all_lengths(:, i))]);
end
% Calculate length using the median of max lengths in all trials
lengths = median(al_sub_min);

%% Pre-processing

for k = 1:length(files) % Each trial
    cd(data_folder);
    rtf = dir(['Trial-', num2str(k), '.mat']);
    load(rtf.name);
    split = strsplit(rtf.name, '.');
    DATAID = split{1};
    count = 0;
    
    % Remove tags
    for j = new_tag_ind
        count = count + 1; 
        XPlot_new(count, :) = XPlot(j, :);
        YPlot_new(count, :) = YPlot(j, :);
        ZPlot_new(count, :) = ZPlot(j, :);
        RollPlot_new(count, :) = RollPlot(j, :);
        PitchPlot_new(count, :) = PitchPlot(j, :);
        YawPlot_new(count, :) = YawPlot(j, :);
    end
    XPlot = XPlot_new;
    YPlot = YPlot_new;
    ZPlot = ZPlot_new;
    RollPlot = RollPlot_new;
    PitchPlot = PitchPlot_new;
    YawPlot = YawPlot_new;

    if ~isempty(st_ed)
        sted = st_ed(k, :);
    else
        sted = [];
    end

    %% ========Correction and smoothing=======
    num_tags = size(XPlot, 1);
    seglets = num_tags - 1;

    % Correct orientation
    for i = 1:num_tags % Each tag
        for j = 1:length(XPlot) % Each frame
            R = EULERZYX_int([RollPlot(i,j), PitchPlot(i,j), YawPlot(i,j)]);
            R_corrected = R*Rcrt{i};
            crtd_ang = EULERZYXINV_int(R_corrected);
            RollPlot(i,j) = crtd_ang(1);
            PitchPlot(i,j) = crtd_ang(2);
            YawPlot(i,j) = crtd_ang(3);
        end
    end

    for i = 1:num_tags % each tag
        % Containers for recording axis unit vectors
        e1s = zeros(3, length(XPlot));
        e2s = zeros(3, length(XPlot));
        e3s = zeros(3, length(XPlot));
        % Record axis unit vectors
        for j = 1:length(XPlot) % each frame
            R = EULERZYX_int([RollPlot(i,j), PitchPlot(i,j), YawPlot(i,j)]);
            e1s(1:3, j) = R(1:3, 1);
            e2s(1:3, j) = R(1:3, 2);
            e3s(1:3, j) = R(1:3, 3);
        end
        % Smooth each unit vector temporally, this does not close in SO(3)
        for k = 1:3 
            e1s(k,:) = smooth2a(e1s(k,:),0,4);
            e2s(k,:) = smooth2a(e2s(k,:),0,4);
            e3s(k,:) = smooth2a(e3s(k,:),0,4);
        end
        % Revert to euler angles for recording
        for j = 1:length(XPlot) % Each frame
            tempR = [e1s(1:3, j), e2s(1:3, j), e3s(1:3,j)];
            if isnan(tempR(1,1))
                RollPlot(i,j) = nan;
                PitchPlot(i,j) = nan;
                YawPlot(i,j) = nan;
            else
                [U,~,V] = svd(tempR);
                ang = EULERZYXINV_int(U*V');
                RollPlot(i,j) = ang(1);
                PitchPlot(i,j) = ang(2);
                YawPlot(i,j) = ang(3);
            end
        end
    end

    %% Apply tapering translation to bring tags down to body
    tot_lgth = sum(lengths);
    tl = cumsum(lengths)/tot_lgth; % Ratio milestone for each tag (2:end)
    tli = [1, round(100*tl)]; % Convert to percent
    tp_r = taper(tli);

    for i = 1:num_tags
        for j = 1:length(XPlot)
            R = EULERZYX_int([RollPlot(i, j), PitchPlot(i, j), YawPlot(i, j)]);
            trans = [XPlot(i,j); YPlot(i,j); ZPlot(i,j)];
            trans = trans - (th+ tp_r(i)*rad)*R(1:3, 3);
            XPlot(i, j) = trans(1);
            YPlot(i, j) = trans(2);
            ZPlot(i, j) = trans(3);
        end
    end


    %% Save Data, you can modify the save name here
    cd(save_folder)
    save([DATAID, '_pp_kt-', kill_str, '.mat'], 'Rcrt', 'lengths', 'XPlot', 'YPlot', 'ZPlot', 'RollPlot', 'PitchPlot', 'YawPlot', 'sted', 'order');

    clear XPlot YPlot ZPlot RollPlot PitchPlot YawPlot XPlot_new YPlot_new ZPlot_new RollPlot_new PitchPlot_new YawPlot_new 
end