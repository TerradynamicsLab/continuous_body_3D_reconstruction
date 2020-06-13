%% Jin Seob Kim 2006, Thomas W.Mitchel 2018 tmitchel@jhu.edu
% Modified by Qiyuan:
%   added comments and fixed bugs; can limit frame range (e.g. for a quick validation)
%% Reconstruct the snake body from BEEtag data
% Inputs:

%     tag_data: File name of pre-processed tag data you want to reconstruct
%
%     save_name: Character string which defines what you'll save the data as.
%     The data will be saved with the name [save_name, '.mat'];
%
%     RML: A 1x3 array consisting of the most basic measurments of the snake
%     RML(1) = radius of snake at midpoint in mm
%     RML(2) = mass of snake in grams
%     RML(3) = total length of snake from head to tail in mm

%     frameRange: [start, ending] frame to reconstruct, [] if for all the
%     frames

%     lgth_adj_p: a vector, each element is a proportion to adjust length

% Optional Arguements:

%   E: A 2x1 vector containing the Youngs moduli corresponding to the bending
%   stiffness
%       E(1) = E1, the Young's modulus corresponding to side to side bending
%       E(2) = E2, the Young's modulus corresponding to up and down bending

%   G: A scalar which defines the shear modulus

% If E and G are not defined, the reconstruction will assume that that the
% snake is isoparameteric and use moduli corresponding to fish and human
% tissue


function reconstruct_body_elastic(tag_data, save_name, RML, frameRange, lgth_adj_p, E, G)
%% Variable definition/intialization/preparation
global N L a10 a20 c0 s1 s2 d0 m_s_g g_Iap K Kinv kv epsilon  B

ER_4RK = @extrod_4RK_F; % Return the distal frame given eta(0)
ER_4RK_v1 = @extrod_4RK_F_all; % Return all frames and twists of elements on the segment given eta(0)

%% Thresholds

% These two detect when a conformation is not 'snakelike', I would strongly
% reccomend using default values
writhe_thresh_max = .4; % Default = .4
twist_thresh_max = .4; % Default = .4

% Congergence threshold, iterative solution stops when error is below this
% number, weighted Frank Park norm dimensionless units.
conv_thresh = .1;

%% Arclength point number
% How many elements you want per segment, I'd reccomend this to be at least
% 700, any lower will lead to worse-fitting backbone

N = 700;

%% Load marker tracking data

load(tag_data)

%% ======= Arc length and time count initialization =======%

% 8 initial guesses, each row is a xi with bending in different directions and extension
N0 = [[.01, 0, 0;...
    -.01, 0, 0;...
    0, .01, 0;...
    0, -.01, 0;...
    .01, .01, 0;...
    -.01, .01, 0;...
    .01, -.01, 0;...
    -.01, -.01, 0], zeros(8, 2), ones(8, 1)];
seglets = size(XPlot,1)-1; % (Qiyuan) assume seglegts is the number of section
ig = zeros(6, seglets); % Eta, (undetermined coefficients \mu and \lambda)
for j = 1:seglets % seglets, what is it?
    ig(:, j) = N0(1, :)'; % Initialize with the same value (1st initial guess)
end

% For time evaluation
total_time_seg = 0;
sg_cnt = 0;

%%
%=============================%
%%%%%%%% Physical data %%%%%%%%
%=============================%

%% Radius, mass, length

% Radius of the snake in mm
rad = RML(1);

% Mass of the snake in grams
s_m = RML(2);

% Total body length of the snake in mm
tot_lgth = RML(3);

%% Compute some important quantities
s_m = .001*s_m; %convert to kg
rad_m = .001*rad; %radius of snake in m
tot_lgth = .001*tot_lgth; % length of snake in m

s_A = pi*(rad_m)^2; %cross sectional area of snake in m^2
s_V = tot_lgth*s_A; %total volume of snake m^3
den = s_m/s_V; %density of snake body in kg/m^3
grav = 9.8; %gravitational acceleration in m/s^2

m_s_g = -s_A*den*grav*.001;

%% Define material properties

% If no specific moduli are given, default to isotropic case with fish +
% human tissue parameters

if nargin < 6
    E = (10^5); %Youngs modulus for tissue in N/m^2, ranges from 10^5 ~ 10^6
    E1 = E*(10^-6); %Convert to N/mm^2
    E2 = E1;
    nu = .29;%.47; %Poisson's ratio for soft tissue, .45 < nu <.49
    G = E1/(2*(1 + nu)); %modulus of rigidity in N/mm^2
else
    if length(E) == 1
        E1 = E(1);
        E2 = E(2);
    else
        E1 = E;
        E2 = E;
    end
    
    %Convert all to N/mm^2
    E1 = (1e-6)*E1;
    E2 = (1e-6)*E2;
    G = (1e-6)*G;
end

%% Define stiffness coefficents

I = (pi/4)*(rad^4); %area  moment of intertia
J = (pi/2)*(rad^4); %Torsion constant
cs_A = pi*(rad^2); %Cross section

a10 = E1*I; %side to side bending stiffness
a20 = E2*I;
c0 = G*J; %Torsion stiffness
d0 = mean([E1, E2])*cs_A;
s1 = G*cs_A;
s2 = G*cs_A;

%% Define internal stiffness matrix

B=[a10,0,  0;...
    0,  a20,0;...
    0,  0,  c0];
C=[0,0,0;...
    0,0,0;...
    0,0,0];
D=[s1,0,0;...
    0,s2,0;...
    0,0,d0];
K=[B,C;C',D];
kv=K*[zeros(3, 1);0;0;1];

if cond(K) > 1e10
    Kinv = pinv(K, 0.0000001);
else
    Kinv=inv(K);
end

count = 0;
segm = 0;


%% =================== Begin frame by frame reconstruction =====================
if isempty(frameRange)
    frameRange = [1, size(XPlot,2)];
end

% Initialize container
% Segments variable (record using g in SE(3))
segments = repmat(struct('backbone',{cell(size(XPlot,2),1)},...
    'bodyvel',{cell(size(XPlot,2),1)}),seglets);
% Seg variable (record using X/Y/Z/Yaw/Pitch/Roll instead of g)
seg = repmat(struct('B_L',nan,'x',[],'y',[],'z',[],'roll',[],'pitch',[],'yaw',[]),...
    seglets);
for indFrame = frameRange(1):frameRange(2)
    count = count + 1;
    ind = indFrame;
    for indSeg = 1:seglets
        % Length of the segment
        L = lengths(indSeg)*(1+lgth_adj_p(indSeg));
        seg(indSeg).B_L = L;
        sg_cnt = sg_cnt + 1;
        tic
        disp(['>>> frame # =', num2str(indFrame), ', seg # =', num2str(indSeg)])
        if ~isnan(XPlot(indSeg, indFrame)) && ~isnan(XPlot(indSeg+1, indFrame))
            segm = segm + 1;
            tag1 = indSeg;
            tag2 = indSeg + 1;
            
            %% Define end constraints imposed by tags
            
            g_Ia = [EULERZYX_int([RollPlot(tag1, ind), PitchPlot(tag1, ind), YawPlot(tag1, ind)]), [XPlot(tag1, ind); YPlot(tag1, ind); ZPlot(tag1, ind)]; zeros(1,3), 1];
            g_Ib = [EULERZYX_int([RollPlot(tag2, ind), PitchPlot(tag2, ind), YawPlot(tag2, ind)]), [XPlot(tag2, ind); YPlot(tag2, ind); ZPlot(tag2, ind)]; zeros(1,3), 1];
            
            %% Transform into algorithm coordinates
            g_Iap = [expm(hat(g_Ia(1:3, 2))*(-pi/2))*g_Ia(1:3, 1:3), g_Ia(1:3, 4); zeros(1,3), 1];
            g_Ibp = [expm(hat(g_Ib(1:3, 2))*(-pi/2))*g_Ib(1:3, 1:3), g_Ib(1:3, 4); zeros(1,3), 1];
            
            gd = inv_SE(g_Iap)*g_Ibp; % transform from a to b
            
            gdOG = gd;
            
            ds = L/N; % Length of each small arc element
            
            %% SERIOUSLY DON"T MODIFY ANYTHING BELOW HERE UNLESS
            
            %% =========Initialize Tries========
            
            if count == 1 % First frame
                t_start = 3;
            else
                t_start = 1;
            end
            
            fnd = 0; % found solution flag
            skip = 0;
            n_elst = NaN*zeros(6, 1 + 8);
            
            %% Begin Tries
            for tries = t_start:3 % First frame, try once; otherwise 3 times
                if tries < 3
                    if tries == 1
                        M = 50; % Number of sections in the artificial end pose trajectory / Number of steps
                        epsilon = 1e-7;
                        dt = 1/M;
                        n0 = ig(1:6, indSeg);
                    else % tries == 2
                        M = 210;
                        epsilon = 1e-10;
                        dt = 1/2000;
                        n0 = ig(1:6, indSeg);
                    end
                    
                    if skip == 0
                        
                        gs=zeros(4*(N+1),4*(M+1)); % collection of snake conformation
                        g_dist=zeros(4*(M+1),4); % collection of distal end poses
                        
                        %% ===== calculating the 1st conformation
                        j=1;
                        
                        [gt0,~]=ER_4RK_v1(n0, L, N, K, Kinv, kv); % the 1st conformation (g of all the small elegments in this seglet)
                        
                        g0=gt0(4*N+1:4*(N+1),:); % distal end of the 1st conformation
                        g_dist(4*(j-1)+1:4*j,:)=g0;
                        gs(:,4*(j-1)+1:4*j)=gt0;
                        
                        %% ======= information for defining gp =======%
                        Ad=gd(1:3,1:3);
                        avd=gd(1:3,4);
                        
                        %% ======= variational calculus and inverse kinematics =======%
                        
                        %% UNDER NO CIRCUMSTANCES MODIFY THE BELOW
                        
                        n=zeros(6,M+1); % Eta for each step
                        
                        %=======================%
                        %%%%%%% iteration %%%%%%%
                        %=======================%
                        g=g0; n(:,1)=n0; % initialize
                        met_path_check = zeros(1, M+1);
                        
                        for j=2:M+1 % At each step, push the distal frame a bit to the desired one
                            t=(j-1)*dt;
                            
                            %======= Jacobian =======%
                            %%% calculating Jacobian on each direction
                            dg1=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[1;0;0;0;0;0], L, N, K, Kinv, kv)-g);
                            dg2=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;1;0;0;0;0], L, N, K, Kinv, kv)-g);
                            dg3=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;0;1;0;0;0], L, N, K, Kinv, kv)-g);
                            dg4=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;0;0;1;0;0], L, N, K, Kinv, kv)-g);
                            dg5=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;0;0;0;1;0], L, N, K, Kinv, kv)-g);
                            dg6=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;0;0;0;0;1], L, N, K, Kinv, kv)-g);
                            
                            g_verse = inv_SE(g);
                            H1 = g_verse*dg1; H2 = g_verse*dg2; H3 = g_verse*dg3;
                            H4 = g_verse*dg4; H5 = g_verse*dg5; H6 = g_verse*dg6;
                            
                            H1=[0.5*(H1(1:3,1:3)-H1(1:3,1:3)'),H1(1:3,4);zeros(1,3),0];
                            H2=[0.5*(H2(1:3,1:3)-H2(1:3,1:3)'),H2(1:3,4);zeros(1,3),0];
                            H3=[0.5*(H3(1:3,1:3)-H3(1:3,1:3)'),H3(1:3,4);zeros(1,3),0];
                            H4=[0.5*(H4(1:3,1:3)-H4(1:3,1:3)'),H4(1:3,4);zeros(1,3),0];
                            H5=[0.5*(H5(1:3,1:3)-H5(1:3,1:3)'),H5(1:3,4);zeros(1,3),0];
                            H6=[0.5*(H6(1:3,1:3)-H6(1:3,1:3)'),H6(1:3,4);zeros(1,3),0];
                            
                            for k=1:6 % Calculate "right" dual vector of each Hk
                                eval(['Jr(:,' num2str(k) ')=vect(H' num2str(k) ');'])
                                % Seems that this eval is just to prevent error (could be replaced by try catch after modifying Hk data structure)
                            end
                            
                            %Check that Jr exists lol
                            
                            uhoh = 0;
                            for i = 1:6
                                if isnan(Jr(1,i))
                                    uhoh = uhoh +1;
                                end
                            end
                            if uhoh ~=0 % Some Jr doesn't exist
                                print('ERROR: invalid tag data. Check that the rotational components of the input data are actually in SO(3)');
                                print('If this error keeps reoccuring, please contact Tommy at tmitchel@jhu.edu');
                                return
                            end
                            
                            
                            %======= inverse problem solving procedure =======%
                            %%% artificial path of the distal frame
                            A=g(1:3,1:3);
                            av=g(1:3,4);
                            gp=[A*exp_so((t-dt)*log_SO(A'*Ad)),(t-dt)*(avd-av)+av;zeros(1,3),1];
                            gpprime=[A*log_SO(A'*Ad)*exp_so((t-dt)*log_SO(A'*Ad)),avd-av;zeros(1,4)];
                            
                            %%% velocity track, Eqn. (14)
                            if cond(Jr) > 1e10
                                dEta=pinv(Jr,0.0000001)*Adjoint(inv_SE(g)*gp)*vect(inv_SE(gp)*gpprime);
                            else
                                dEta=Jr\Adjoint(inv_SE(g)*gp)*vect(inv_SE(gp)*gpprime);
                            end
                            n1=n(:,j-1)+dt*dEta;
                            
                            %%% position track (correction term), Eqn. (15)
                            
                            g_cor = inv_SE(g)*gp;
                            X_cor=log_SE(g_cor);
                            if cond(Jr) > 1e10
                                n2=pinv(Jr,0.0000001)*Adjoint(inv_SE(g)*gp)*vect(X_cor);
                            else
                                n2=Jr\Adjoint(inv_SE(g)*gp)*vect(X_cor);
                            end
                            
                            n(:,j)=real(n1+n2); % Update eta for the j-th step
                            
                            %%% actual conformation for j-th step
                            [gt,~]=ER_4RK_v1(n(:,j), L, N, K, Kinv, kv);
                            g=gt(4*N+1:4*(N+1),:); % distal end value of each conformation
                            g_dist(4*(j-1)+1:4*j,:)=g;
                            gs(:,4*(j-1)+1:4*j)=gt;
                            
                            % for display
                            gp2=[A*exp_so((t)*log_SO(A'*Ad)),(t)*(avd-av)+av;zeros(1,3),1];
                            metric_path=sqrt(norm(g(1:3,1:3)-gp2(1:3,1:3),2)^2+norm(g(1:3,4)-gp2(1:3,4))^2);
                            met_path_check(j) = metric_path;
                            
                            % Quit trying conformation if unlikely to converge
                            if j > 65 % Over 65 steps
                                if min(met_path_check(1, j-6:j)) > 1.75|| abs(met_path_check(1, j) - met_path_check(1, j-1)) > 8 || met_path_check(j) > 8
                                    % And it's still far away from the desired distal frame
                                    metric = 100;
                                    break
                                end
                            end
                            
                            % distance from g to gd
                            metric=sqrt(norm(g(1:3,1:3)-gd(1:3,1:3),2)^2+norm(g(1:3,4)-gd(1:3,4))^2);
                            if isnan(metric)
                                metric = 100;
                                break
                            end
                            if metric <= conv_thresh
                                Ms = j;
                                break
                            end
                            
                        end
                        
                        
                        %% Evaluate conformation to check if it is snake-like
                        
                        if metric <= conv_thresh % Distal frame close enough to the desired one
                            
                            [gcan,xican] = ER_4RK_v1(n(:, Ms), L, N, K, Kinv, kv); % Calculate seglet deformation now
                            rythe = writhe(gcan, L);
                            wcan = xican(1:3, :);
                            wist = twist_n(wcan, L);
                            
                            if abs(wist) < twist_thresh_max && abs(rythe) <= writhe_thresh_max % No significant writhe or twist
                                n_elst(:, 1) =  n(:, Ms);
                                ig(1:6, indSeg) = n(:, Ms); % Record eta (undetermined coefficients) of this seglet
                                gF = gcan; % Record conformation
                                wF = wcan; % Record w of each small element
                                fnd = 1;
                                break
                            else
                                if tries == 1
                                    skip = 1;
                                    metric = 100;
                                end
                            end
                        end
                        
                    end % This end should be the one that ends the if skip = 0 statement
                    
                    %% If previous solutions fail to converge to conformation, then try the 8 initial guesses
                else % Tries == 3
                    
                    M = 200; % Max iteration number (number of steps)
                    epsilon = 1e-10;
                    dt = 1/2000;
                    gd = gdOG;
                    nf = zeros(8,6);
                    
                    %% Reorder intitial guesses based on distance from metric
                    IGmetric = zeros(8, 1);
                    % Frank-Parker distance between the endpoint reconstructed and the next tag
                    for qq = 1:8 % For each initial guess
                        [g0_CHK, ~]=ER_4RK_v1(N0(qq, :)', L, N, K, Kinv, kv); % g of each element on the reconstructed rod
                        IGmetric(qq) =  sqrt(norm(g0_CHK(4*N+1:4*N+3,1:3)-gd(1:3,1:3),2)^2+norm(g0_CHK(4*N+1:4*N+3,4)-gd(1:3,4))^2);
                    end
                    % Sort distance from small to large
                    IG_reorder = zeros(8, 1); % Index of 8 guesses, small distance to large
                    for qq = 1:8
                        [~, indCHK] = min(IGmetric(:, 1));
                        IG_reorder(qq) = indCHK;
                        IGmetric(indCHK) = NaN;
                    end
                    
                    
                    %% Try initial gueses in order of most likley to succeed in producing a valid conformation
                    % The calculation is basically the same as previous
                    % tries
                    %% DON'T MODIFY ANYTHING BELOW HERE
                    for indGuess = 1:8
                        met_path_check = zeros(1, M+1);
                        metric = 100;
                        n0 = N0(IG_reorder(indGuess), :)';
                        
                        gs=zeros(4*(N+1),4*(M+1)); % collection of DNA conformation
                        g_dist=zeros(4*(M+1),4); % collection of distal end poses
                        
                        %% ===== calculating the 1st conformation
                        j=1; % Iteration index
                        [gt0,~]=ER_4RK_v1(n0, L, N, K, Kinv, kv); % the 1st conformation
                        
                        g0=gt0(4*N+1:4*(N+1),:); % distal end of the 1st conformation
                        g_dist(4*(j-1)+1:4*j,:)=g0;
                        gs(:,4*(j-1)+1:4*j)=gt0;
                        
                        %======= information for defining gp =======%
                        Ad=gd(1:3,1:3); % Desired rotation from this to next tag
                        avd=gd(1:3,4); % Desired translation from this to next tag
                        
                        %% ======= variational calculus and inverse kinematics =======%
                        n=zeros(6,M+1);
                        n(:,1)=n0;
                        
                        %=======================%
                        %%%%%%% iteration %%%%%%%
                        %=======================%
                        g=g0; n(:,1)=n0; % initialize
                        
                        for j=2:M+1 % Iterate
                            t=(j-1)*dt;
                            
                            %======= Jacobian =======%
                            %%% calculating Jacobian
                            dg1=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[1;0;0;0;0;0], L, N, K, Kinv, kv)-g);
                            dg2=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;1;0;0;0;0], L, N, K, Kinv, kv)-g);
                            dg3=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;0;1;0;0;0], L, N, K, Kinv, kv)-g);
                            dg4=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;0;0;1;0;0], L, N, K, Kinv, kv)-g);
                            dg5=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;0;0;0;1;0], L, N, K, Kinv, kv)-g);
                            dg6=1/epsilon*(ER_4RK(n(:,j-1)+epsilon*[0;0;0;0;0;1], L, N, K, Kinv, kv)-g);
                            
                            g_verse = inv_SE(g);
                            H1 = g_verse*dg1; H2 = g_verse*dg2; H3 = g_verse*dg3;
                            H4 = g_verse*dg4; H5 = g_verse*dg5; H6 = g_verse*dg6;
                            
                            H1=[0.5*(H1(1:3,1:3)-H1(1:3,1:3)'),H1(1:3,4);zeros(1,3),0];
                            H2=[0.5*(H2(1:3,1:3)-H2(1:3,1:3)'),H2(1:3,4);zeros(1,3),0];
                            H3=[0.5*(H3(1:3,1:3)-H3(1:3,1:3)'),H3(1:3,4);zeros(1,3),0];
                            H4=[0.5*(H4(1:3,1:3)-H4(1:3,1:3)'),H4(1:3,4);zeros(1,3),0];
                            H5=[0.5*(H5(1:3,1:3)-H5(1:3,1:3)'),H5(1:3,4);zeros(1,3),0];
                            H6=[0.5*(H6(1:3,1:3)-H6(1:3,1:3)'),H6(1:3,4);zeros(1,3),0];
                            for k=1:6 % Each component/direction
                                eval(['Jr(:,' num2str(k) ')=vect(H' num2str(k) ');'])% Right Jacobian
                            end
                            
                            % Check that Jr exists
                            uhoh = 0;
                            for i = 1:6
                                if isnan(Jr(1,i))
                                    uhoh = uhoh +1;
                                end
                            end
                            if uhoh ~=0
                                print('ERROR: invalid tag data. Check that the rotational components of the input data are actually in SO(3)');
                                print('If this error keeps reoccuring, please contact Tommy at tmitchel@jhu.edu');
                                return
                            end
                            
                            
                            %======= inverse problem solving procedure =======%
                            % Refer to Conformational analysis of stiff
                            % chiral polymers with end-constraints for
                            % detailed explanations and figrues
                            
                            %%% artificial path from current distal frame
                            %%% to the desired distal frame
                            A=g(1:3,1:3);
                            av=g(1:3,4);
                            gp=[A*exp_so((t-dt)*log_SO(A'*Ad)),(t-dt)*(avd-av)+av;zeros(1,3),1];
                            gpprime=[A*log_SO(A'*Ad)*exp_so((t-dt)*log_SO(A'*Ad)),avd-av;zeros(1,4)];
                            
                            %%% velocity track, Eqn. (14)
                            if cond(Jr) > 1e10
                                dEta=pinv(Jr,0.0000001)*Adjoint(inv_SE(g)*gp)*vect(inv_SE(gp)*gpprime);
                            else
                                dEta = Jr\Adjoint(inv_SE(g)*gp)*vect(inv_SE(gp)*gpprime);
                            end
                            n1=n(:,j-1)+dt*dEta;
                            
                            %%% position track (correction term), Eqn. (15)
                            
                            g_cor = inv_SE(g)*gp;
                            X_cor=log_SE(g_cor);
                            if cond(Jr) > 1e10
                                n2=pinv(Jr,0.0000001)*Adjoint(inv_SE(g)*gp)*vect(X_cor);
                            else
                                n2 = Jr\Adjoint(inv_SE(g)*gp)*vect(X_cor);
                            end
                            
                            n(:,j)=real(n1+n2);
                            
                            %%% actual conformation for j-th step
                            [gt,~]=ER_4RK_v1(n(:,j), L, N, K, Kinv, kv);
                            g=gt(4*N+1:4*(N+1),:); % distal end value of each conformation
                            g_dist(4*(j-1)+1:4*j,:)=g;
                            gs(:,4*(j-1)+1:4*j)=gt;
                            
                            % for display
                            gp2=[A*exp_so((t)*log_SO(A'*Ad)),(t)*(avd-av)+av;zeros(1,3),1];
                            metric_path=sqrt(norm(g(1:3,1:3)-gp2(1:3,1:3),2)^2+norm(g(1:3,4)-gp2(1:3,4))^2);
                            met_path_check(j) = metric_path;
                            
                            % Abort guess if it appears to be diverging
                            if j > 65
                                if min(met_path_check(1, j-6:j)) > 1.75|| abs(met_path_check(1, j) - met_path_check(1, j-1)) > 8 || met_path_check(j) > 8
                                    metric = 100;
                                    break
                                end
                            end
                            
                            % Check for convergence
                            metric = sqrt(norm(g(1:3,1:3)-gd(1:3,1:3),2)^2+norm(g(1:3,4)-gd(1:3,4))^2);
                            if isnan(metric)
                                metric = 100;
                                break
                            end
                            if metric <= conv_thresh
                                nf(indGuess, 1:6) = n(:, j)';
                                break
                            end
                            
                        end
                        
                        %% Evaluate conformation if converged
                        if metric > conv_thresh % After all these trials no success
                            nf(indGuess, 1:6) = [NaN, NaN, NaN, NaN, NaN, NaN];
                        else
                            [gtIG, xiIG] = ER_4RK_v1(nf(indGuess, 1:6)', L, N, K, Kinv, kv);
                            wIG = xiIG(1:3, :);
                            rythe = writhe(gtIG, L);% Qunatifies the amount of coiling
                            wist = twist_n(wIG, L);
                            
                            if abs(wist) < twist_thresh_max && abs(rythe) <= writhe_thresh_max
                                n_elst(:, (tries-3)*8 + indGuess + 1) =  nf(indGuess, 1:6)';
                                ig(1:6, indSeg) = nf(indGuess, 1:6)'; % Update the initial guess with the current one
                                gF = gtIG;
                                wF = wIG;
                                fnd = 1;
                                break
                            end
                            
                        end
                    end
                    
                    %% If you've found a snake-like conformation, exit the solution process
                    if fnd ~= 0
                        break
                    end
                    
                    
                end % end if tries < 3
                
            end %%%%%%%%%% end tries and solution loop
            
            %% Convert found conformation from algorithm coordinates back to snake coordinates
            if fnd ~= 0
                gbb = zeros(4*(N+1), 4);
                
                for j = 1:4*(N+1)
                    if mod(j, 4) == 1
                        gbb(j:j+3, 1:4) = g_Iap*gF(j:j+3, 1:4);
                        gbb(j:j+2, 1:3) = exp_so(hat(gbb(j:j+2, 2))*(pi/2))*gbb(j:j+2, 1:3);
                    end
                end
                
                segments(indSeg).backbone{indFrame} = gbb;
                segments(indSeg).bodyvel{indFrame} = wF;
            end
        end % end of if ~isnan statment
        %% Update time counter
        a = toc;
        total_time_seg = total_time_seg + a;
        avg_t_seg = total_time_seg/sg_cnt;
        total_hrs = total_time_seg/3600;
        disp(['Avg time per segment = ', num2str(avg_t_seg), ' sec'])
        disp(['Time elapsed = ' num2str(total_hrs), ' hrs'])
        
    end %==========end iii = 1:seglets
end %================end ii = 1:step:length(XPlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% Convert segments to seg

%===fill===
% Empty ones are filled with nan

for j = 1:seglets
    f0(j) = NaN;
    ff(j) = NaN;
    for i = 1:size(XPlot, 2)
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


%===Convert===
for i = 1:size(XPlot, 2)
    for j = 1:seglets
        count = 0;
        gseg = segments(j).backbone{i};
        Nseg=N;%Qiyuan assumed this
        for k = 1:4:4*(Nseg + 1)
            count = count + 1;
            if ~isempty(segments(j).backbone{i})
                R = gseg(k:k+2, 1:3);
                ang = EULERZYXINV_int(R);
                seg(j).x(count, i) = gseg(k,4);
                seg(j).y(count, i) = gseg(k+1,4);
                seg(j).z(count, i) = gseg(k+2,4);
                seg(j).roll(count, i) = ang(1);
                seg(j).pitch(count, i) = ang(2);
                seg(j).yaw(count, i) = ang(3);
            else
                seg(j).x(count, i) = NaN;
                seg(j).y(count, i) = NaN;
                seg(j).z(count, i) = NaN;
                seg(j).roll(count, i) = NaN;
                seg(j).pitch(count, i) = NaN;
                seg(j).yaw(count, i) = NaN;
            end
        end
    end
end

st_frame = max(f0);
ed_frame = max(ff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555
%% Save data
% cd(savepath)
save(save_name, 'segments', 'seg', 'st_frame', 'ed_frame','lgth_adj_p', '-v7.3');

end