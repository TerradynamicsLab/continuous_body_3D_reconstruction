% Post-process reconstructed data produced by Body_reconst_elastic.m

% Output data contains 4 variables:
%   st_frame: starting video frame
%   ed_frame: ending video frame
%   seg: n x 1 structure, where n is the number of segments
%        Each element contains 7 fields: B_L is the length of this segment,
%        x/y/z/roll/pitch/yaw are x/y/z coordinates and roll/pitch/yaw
%        Euler angles (ZYX convention)
%   segments: n x 1 structure, where n is the number of segments
%             Each element contains 2 fields: backbone and bodyvel
%             backbone{i}((4*j-3):4*j,:) is the spatial configuration
%             matrix of the j-th finite element in the i-th video frame
%             bodyvel{i}(:,j) is the unit length rotational deformation
%             of the j-th finite element in the i-th video frame

close all;clearvars;clc

%% Input
% Raw reconstructed backbone data path
rawdp = ['D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Data\',...
    'Reconstructed\']; 
% Pre-processed data path
ppdp = ['D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Data\',...
    'Pre-processed\'];  
% Folder to save post-processed data
savepath = ['D:\Dropbox (Terradynamics Lab)\All\Snake Modeling\Codes\To release with sample data\Data\',...
    'Post-processed\'];

% Number of elements in each segment
Nseg = 700; 
N = Nseg;

%% Preparation
cd(rawdp);
files = dir('*_backbone*.mat');

cd(ppdp)
disc = dir('*_pp*.mat');

%% Processing
for FILES = 1:length(files) % each file
    % Load data
    cd(rawdp)
    load(files(FILES).name);
    cd(ppdp)
    load(disc(FILES).name);
    
    seglets = length(segments); % Number of segments in the body
    
    % End frame
    for mu = 1:length(segments)
        efc(mu) = length(segments(mu).backbone);
    end
    endframe = min(efc);    

    % Find out the 1st/last frame that tracked the segment
    for j = 1:seglets
        f0(j) = NaN; % 1st frame
        ff(j) = NaN; % last frame
        for i = 1:endframe
            if isnan(f0(j)) 
                if ~isempty(segments(j).backbone{i})
                    f0(j) = i;
                end
            end
            if isnan(ff(j))
               if ~isempty(segments(j).backbone{length(segments(j).backbone)+ 1 - i})
                   ff(j) = length(segments(j).backbone)+ 1 - i;
               end
            end
            if ~isnan(f0(j)) && ~isnan(ff(j))
                break
            end
        end
    end

    % Concatenate all small segments
    r(1).x = NaN(seglets*(N+1), max(ff));
    r(1).y = NaN(seglets*(N+1), max(ff));
    r(1).z = NaN(seglets*(N+1), max(ff));
    r(1).roll = NaN(seglets*(N+1), max(ff));
    r(1).pitch = NaN(seglets*(N+1), max(ff));
    r(1).yaw = NaN(seglets*(N+1), max(ff));
    for j = 1:seglets
        r(1).x((j-1)*(N+1) + 1: j*(N+1), max(f0):min(ff)) = seg(j).x(1:N+1, max(f0):min(ff));
        r(1).y((j-1)*(N+1) + 1: j*(N+1), max(f0):min(ff)) = seg(j).y(1:N+1, max(f0):min(ff));
        r(1).z((j-1)*(N+1) + 1: j*(N+1), max(f0):min(ff)) = seg(j).z(1:N+1, max(f0):min(ff));
        r(1).roll((j-1)*(N+1) + 1: j*(N+1), max(f0):min(ff)) = seg(j).roll(1:N+1, max(f0):min(ff));
        r(1).pitch((j-1)*(N+1) + 1: j*(N+1), max(f0):min(ff)) = seg(j).pitch(1:N+1, max(f0):min(ff));
        r(1).yaw((j-1)*(N+1) + 1: j*(N+1), max(f0):min(ff)) = seg(j).yaw(1:N+1, max(f0):min(ff));
    end
    
     % === fill missing temporally ===   
    for u = 1:seglets*(Nseg + 1)
        g2fill = zeros(4, 4, length(max(f0):min(ff)));
        % Convert to g
        for i = max(f0):min(ff) % Frame
            g2fill(:, :, i) = [EULERZYX_int([r(1).roll(u, i), r(1).pitch(u, i), r(1).yaw(u, i)]), [r(1).x(u, i); r(1).y(u, i); r(1).z(u, i)]; zeros(1, 3), 1];
        end
        % Temporally fill missing
        g2fill = fillmissing_SE_3(g2fill);
        % Convert back to x/y/z/roll/pitch/yaw
        for i = max(f0):min(ff)
            r(1).x(u, i) = g2fill(1, 4, i);
            r(1).y(u, i) = g2fill(2, 4, i);
            r(1).z(u, i) = g2fill(3, 4, i);
            ang = EULERZYXINV_int(g2fill(1:3, 1:3, i));
            r(1).roll(u, i) = ang(1);
            r(1).pitch(u, i) = ang(2);
            r(1).yaw(u, i) = ang(3);
        end
    end
    
    % === Simple smoothing of position ===  

    % Position: spatial temporal smoothing
        r(1).x(:, max(f0):min(ff)) = smooth2a(r(1).x(:, max(f0):min(ff)), 2, 5);
        r(1).y(:, max(f0):min(ff)) = smooth2a(r(1).y(:, max(f0):min(ff)), 2, 5);
        r(1).z(:, max(f0):min(ff)) = smooth2a(r(1).z(:, max(f0):min(ff)), 2, 5);
        
    % Convert back to store data
    % seg: x/y/z/pitch/yaw/roll
    for j = 1:seglets  
        seg(j).x(1:N+1, max(f0):min(ff)) = r(1).x((j-1)*(N+1) + 1:j*(N+1), max(f0):min(ff));
        seg(j).y(1:N+1, max(f0):min(ff)) = r(1).y((j-1)*(N+1) + 1:j*(N+1), max(f0):min(ff));
        seg(j).z(1:N+1, max(f0):min(ff)) = r(1).z((j-1)*(N+1) + 1:j*(N+1), max(f0):min(ff));
        seg(j).roll(1:N+1, max(f0):min(ff)) = r(1).roll((j-1)*(N+1) + 1:j*(N+1), max(f0):min(ff));
        seg(j).pitch(1:N+1, max(f0):min(ff)) = r(1).pitch((j-1)*(N+1) + 1:j*(N+1), max(f0):min(ff));
        seg(j).yaw(1:N+1, max(f0):min(ff)) = r(1).yaw((j-1)*(N+1) + 1:j*(N+1), max(f0):min(ff));
    end
    % segments: g
    for j = 1:seglets
        for i = max(f0):min(ff)
            count = 0;
            for k = 1:4:4*(N+1)
                count = count + 1;
                R = EULERZYX_int([seg(j).roll(count, i), seg(j).pitch(count, i), seg(j).yaw(count, i)]);
                segments(j).backbone{i}(k:k+3, 1:4) = [R, [seg(j).x(count, i); seg(j).y(count,i); seg(j).z(count, i)]; zeros(1,3), 1];
            end
        end
    end
       
    st_frame = max(f0);
    ed_frame = max(ff);
    
    % Save data
    cd(savepath)
    splt = strsplit(files(FILES).name, '_backbone');
    save([splt{1}, '_post.mat'],'segments','seg', 'st_frame', 'ed_frame', '-v7.3');
    clear seg segments r g2fill        
end