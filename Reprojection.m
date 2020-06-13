% Visualize segments reconstructed using elastic rod theory and optimization
% QF
close all;clearvars;clc

%% Input
% Folder of the reconstructed data / post-processed data
filepath = ['D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Data',...
    '\Post-processed'];
trialsToPlot = 1:11; % Which trials to re-project

% Folder with the DLT calibration file
calibrationfilepath = 'T:\08-21-17\K 1\2 in\cal02';
videofilepath = 'T:\08-21-17\K 1\2 in';
savePath = 'D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Demo';

%% Reprojection
% Load calibration file
cd(calibrationfilepath)
calFile = dir('cal*.csv');
c = csvread(calFile.name);

for indTrial = trialsToPlot % Trials to plot
    files = dir([filepath,'\Trial-', num2str(indTrial), '_post.mat']);
    if isempty(files)
        continue
    end
    load([files.folder,'\',files.name]);
               
    % Organize data to plot
    XPlot = [];YPlot = [];ZPlot = [];
    for indFrame = 1:length(segments(1).backbone)
        XYZ_this = [];
        for s = 1:length(segments)
            if ~isempty(segments(s).backbone{indFrame})
                XYZ_this = [XYZ_this;...
                    segments(s).backbone{indFrame}(1:4:end,4),...
                    segments(s).backbone{indFrame}(2:4:end,4),...
                    segments(s).backbone{indFrame}(3:4:end,4)];
            else
                XYZ_this = [XYZ_this;...
                    nan(701,3)];
            end
        end
        XPlot = [XPlot,XYZ_this(:,1)];
        YPlot = [YPlot,XYZ_this(:,2)];
        ZPlot = [ZPlot,XYZ_this(:,3)];
    end
    
    % Read raw videos
    vidFile = dir([videofilepath, '\Trial-',num2str(indTrial),'_*\RawVideos\Camera_4_*.avi']);
    vid2 = VideoReader([vidFile.folder, '\', vidFile.name]);    
    
    % Prepare a figure
    close all;
    f1 = figure;
    %set(f1,'position',[10 10 2592 2048]);
    pause(0.01);
    
    % Prepare to write videos
    cd(savePath);
    saveVidName = [files.name(1:end-4),'.avi'];
    vidObj = VideoWriter(saveVidName);vidObj.FrameRate = 100;
    open(vidObj);
    
    % Start writing
    nn=0;start = 1;
    for i = start:size(XPlot,2)
        nn=nn+1;        
        set(0,'currentfigure',f1);hold on;
        
        % Calculate 2-D reprojection points
        uv2Snake = dlt_inverse(c(:,4),[XPlot(:,i),YPlot(:,i),ZPlot(:,i)]);
        
        % Display raw video frame
        I = read(vid2,i);
        if i==start
            I2 = imshow(I);
        else
            I2.CData = I;
        end
        
        % Plot reprojections
        if i>1
            h1.XData = nan;h1.YData = nan;
        end
        hold on;
        h1=plot(uv2Snake(:,1),size(I,1)-uv2Snake(:,2),'color','y','linewidth',4);hold on;
        
        xlim([0 2596]);ylim([0 2048]);

        drawnow
        currFrame = getframe(f1);
        writeVideo(vidObj, currFrame);
        pause(0.05);
        
        clear uv2Snake;
    end
    close(vidObj);
    clear XPlot
    
    set(0,'currentfigure',f1);hold on;
    set(gca,'color','none');axis equal;    
end