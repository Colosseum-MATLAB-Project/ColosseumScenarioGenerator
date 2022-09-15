function [chMatrix_hat, plMatrix_hat] = chApproximation(chMatrix,angleType,scenario)
% Channel approximation function
%
% Author: Miead Tehrani-Moayyed
% Institute for the Wireless Internet of Things, 
% Northeastern University, Boston MA, 02115, USA
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 11-Sep-2022
%
% This function takes chMatrix and apply the ML approximation process
% Input:
%   chMatrix:     channel matrix table 
%   angleType:    can be 'elevation' or 'zenith' angle depending the inputs
%   scenario:     can be 'LOS' or 'NLOS' and defines the delay weight in MCD 
%
% Output:
%   chMatrix_hat    approximated channels of the channel matrix
%                   It is used for Colosseum multi-tap scenario toolchain 
%   plMatrix_hat    path loss matrix of the approximated channels
%                   It is used for Colosseum single-tap scenario toolchain

%% Configuration
% Approximation process parameters
zetaNLOS = 3;                       % MCD distance coefficient
zetaLOS = 6;

% Channel emulator parameters
nTap = 4;                           % Number of non-zero taps
nFIR = 512;                         % Number of FIR filter
fs = 100e6;                         % FIR filter sampling freq.

%% Main process
nTx = size(chMatrix,1);
nRx = size(chMatrix,2);
nSnapshots = size(chMatrix,3);
processData = cell(nTx,nRx,nSnapshots);
chMatrix_hat = cell(nTx,nRx,nSnapshots);
plMatrix_hat = ones(nTx,nRx,nSnapshots) * inf;

for snapshotIdx = 1 : nSnapshots
    for TxIdx = 1 : nTx
        for RxIdx = 1 : nRx

            fprintf('\nSnapshot: %d of %d ,Tx: %d, Rx: %d \n',snapshotIdx,nSnapshots,TxIdx, RxIdx);
            MPC = chMatrix{TxIdx,RxIdx,snapshotIdx}(:,{ 'h' 'tau' 'arrival_theta' 'arrival_phi' 'departure_theta' 'departure_phi'});

            if size(MPC,1) == 0                 %Dealing with no path in the snapshot
                continue
            end

            MPC(MPC.tau>nFIR * 1/fs,:) = [];    %Remove MPC beyond max of excess delay

            if strcmp(angleType,'elevation')    %Adapt the given angle type
                MPC.arrival_theta = 90 - MPC.arrival_theta;
                MPC.departure_theta = 90 - MPC.departure_theta;
            elseif strcmp(angleType,'zenith')
                MPC.arrival_theta = MPC.arrival_theta;
                MPC.departure_theta = MPC.departure_theta;
            else
                error('angleType is invalid. It can only be elevation or zenith');
            end

            if strcmp(scenario,'LOS')
                zeta = zetaLOS;
            elseif strcmp(scenario,'NLOS')
                zeta = zetaNLOS;
            else
                error('Scenario is invalid. It can only be LOS or NLOS');
            end

            % Clustering MPCs
            clusters = MPCclustering(MPC,nTap,zeta);

            % Approximating gains of clustered MPCs
            tau = clusters.estimatedMeans.tau;      
            h = zeros(size(tau));
            gain = zeros(size(tau));
            for idxClass = min(clusters.estimatedLabels):max(clusters.estimatedLabels)
                idxPaths = clusters.estimatedLabels == idxClass;

                h(idxClass) = sum(MPC.h(idxPaths));
                gain(idxClass) = mag2db(abs(h(idxClass)));
                clusters.estimatedMeans.h(idxClass) = h(idxClass);
            end

            clusters.estimatedMeans.h = h;
            clusters.MPC_hat = clusters.estimatedMeans;

            % Re-sampling approximated paths and obtaining FIR Coefficents
            [tau_hat, h_hat] = FIR_coefficient_resampling(tau,h,nFIR,fs);

            ch.clusters = clusters;
            ch.sampled.h = h;
            ch.sampled.tau = tau;
            ch.sampled.h_hat = h_hat;
            ch.sampled.tau_hat = tau_hat;
            ch.sampled.gain = mag2db(abs(sum(h)));
            ch.sampled.gain_hat = mag2db(abs(sum(h_hat)));
            ch.PL = -ch.sampled.gain;
            ch.PL_hat = -ch.sampled.gain_hat;

            chMatrix_hat{TxIdx,RxIdx,snapshotIdx} = table2cell(table(tau_hat,h_hat));
            plMatrix_hat(TxIdx,RxIdx,snapshotIdx) = -ch.sampled.gain_hat;
            processData{TxIdx,RxIdx,snapshotIdx} =  ch;

        end
    end
end

end


%% Functions
function clusters = MPCclustering(MPC,K,zeta)
% Author: Miead Tehrani-Moayyed
% Institute for the Wireless Internet of Things, 
% Northeastern University, Boston MA, 02115, USA
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 11-Sep-2022
%
% Cluster MPCs for one receiver point
% 
% Input
%   MPC - a table of Multi Path Components
%   K   - number of clusters
%   zeta- zeta weight coefficient
%
% Output
%   clusters - output structure
%       MPC - MPC table including the labels column
%       estimatedLabels (Nx1)   - label based on closest mean
%       estimatedMeans (kxd)    - mean estimates
%       vectorMSE (1xnIter)     - MSE as a function of iteration number
%       MSE                     - MSE values
%

params.numberOfClusters = K;
params.numberOfRuns = 10000;    
params.stopTolerance = 1e-8;
params.zeta = zeta;

nMPC = size(MPC,1);
tic
if nMPC > K

    while 1
            [estimatedLabels, estimatedMeans, vectorMSE] = kMeansMCDclustering(MPC, params);

        if FindEmptyCluster( estimatedLabels ) == 0  % if no empty cluster exist
            break
        end
    end
   
    MSE = vectorMSE(end);

else
    estimatedLabels = (1:size(MPC,1))';
    estimatedMeans = MPC;
    vectorMSE = 0;
    MSE = vectorMSE(end);
end
toc

clusters.MPC = MPC;
clusters.MPC.('labels') = estimatedLabels;
clusters.estimatedLabels = estimatedLabels;
clusters.estimatedMeans = estimatedMeans;
clusters.vectorMSE = vectorMSE;
clusters.MSE.values = MSE;

end

function [estimatedLabels, estimatedMeans, vectorMSE] = kMeansMCDclustering(inputData, params)
% Author: Miead Tehrani-Moayyed
% Institute for the Wireless Internet of Things, 
% Northeastern University, Boston MA, 02115, USA
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 11-Sep-2022%
%
% Runs k-means for MPC channel dataset
% Inputs:
%       inputData (Nxd)             - table of MCD data
%       params 
%           .numberOfClusters       - number of clusters
%           .stopTolerance          - abs difference between means for
%                                     stopping criteria
%           .numberOfRuns           - number of random initializations
%           .zeta                   - MCD zeta coefficient
%           .debugMode              - print additional data
%
% Outputs
%       estimatedLabels (Nx1)       - label based on closest mean
%       estimatedMeans (kxd)        - mean estimates
%       vectorMSE (1xnIter)         - MSE as a function of iteration number


if ~isfield(params, 'stopTolerance')
    params.stopTolerance = 1e-8;
end

if ~isfield(params, 'zeta')
    params.zeta=3;
end

if ~isfield(params, 'debugMode')
    params.debugMode=0;
end

numberOfClusters = params.numberOfClusters;
zeta = params.zeta;
stopTolerance = params.stopTolerance;
debugMode = params.debugMode;

% Initialize variables for random runs
runEstimatedLabels = cell(params.numberOfRuns,1);
runMSE = zeros(params.numberOfRuns,1);
runEstimatedMeans = cell(params.numberOfRuns,1);
runVectorMSE = cell(params.numberOfRuns,1);

parfor idxR = 1:params.numberOfRuns
        
    if debugMode==1
        fprintf('On run %d\n', idxR);
    end

    indexRandomInit = randperm(size(inputData,1));
    estimatedLabels = [];

    % Random mean selection
    estimatedMeans = inputData( indexRandomInit(1:numberOfClusters) , :);
    stopCondition = false;
    nIterations = 0;

    while ~stopCondition

        nIterations = nIterations + 1;

        % Initialize
        tmpMeans = estimatedMeans;
        classMSE = zeros(size(inputData, 1), numberOfClusters );

        % Compute MSE for each class
        for idxK = 1:numberOfClusters
            classMSE(:,idxK)=MCD(inputData,estimatedMeans(idxK,:),zeta).^2;
        end

        % Choose the class with smaller MSE
        [MSE, estimatedLabels] = min(classMSE,[],2);

        % Estimate mean for each class
        for idxK = 1:numberOfClusters
            if size(inputData(estimatedLabels==idxK, :),1)>1    % Only Calculate centroids for clusters with multiple members
                estimatedMeans.tau(idxK) = mean(inputData.tau(estimatedLabels==idxK));

                estimatedMeans.arrival_theta(idxK) = meanAngles(inputData.arrival_theta(estimatedLabels==idxK));
                estimatedMeans.arrival_phi(idxK) = meanAngles(inputData.arrival_phi(estimatedLabels==idxK));
                estimatedMeans.departure_theta(idxK) = meanAngles(inputData.departure_theta(estimatedLabels==idxK));
                estimatedMeans.departure_phi(idxK) = meanAngles(inputData.departure_phi(estimatedLabels==idxK));
            end
        end

        % Check the difference with previous estimate
        % diffMeans = min(sum(abs(estimatedMeans-tmpMeans)));
        for idxK = 1:numberOfClusters
            diff = MCD(estimatedMeans,tmpMeans(idxK,:),zeta).^2;
            diffMeans(idxK,1) = diff(idxK);
        end
        diffMeans = mean(diffMeans);

        % Deal with degenerate case of empty class
        if isnan(diffMeans)
            runVectorMSE{idxR}(nIterations) = Inf;
            break
        end

        if mod(nIterations,2) == 1
            if debugMode==1
                fprintf('\t On iteration %d with diffMeans = %e\n', nIterations, diffMeans);
            end
        end

        % Check stopping criteria
        stopCondition = diffMeans < stopTolerance;

        if nIterations>100
            stopCondition=true;
        end

        runVectorMSE{idxR}(nIterations) = mean(MSE);

    end

    runMSE(idxR) = runVectorMSE{idxR}(end);
    runEstimatedLabels{idxR} = estimatedLabels;
    runEstimatedMeans{idxR} = estimatedMeans;

end

% Choose run with smallest MSE
[~, bestRun] = min(runMSE);

estimatedMeans = runEstimatedMeans{bestRun};
estimatedLabels = runEstimatedLabels{bestRun};
vectorMSE = runVectorMSE{bestRun};

end

function d = MCD(MPCi,MPCj,zeta)
% Author: Miead Tehrani-Moayyed
% Institute for the Wireless Internet of Things, 
% Northeastern University, Boston MA, 02115, USA
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 11-Sep-2022
%
% Multipath Component Distance
% Input MPCi,j are tables with the structure below
%   MCD 
%       .tau                -delay
%       .arrival_theta      -AOA zenith
%       .arrival_phi        -AOA azimuth
%       .departure_theta    -AOD zenith
%       .departure_phi      -AOD azimuth
% 
% Output d is MCD distance between MPC i and j
%

if nargin==2
    zeta=3;
end

MCD_AOA=0.5*([sind(MPCi.arrival_theta).*cosd(MPCi.arrival_phi),...
              sind(MPCi.arrival_theta).*sind(MPCi.arrival_phi),...
              cosd(MPCi.arrival_theta)] - ...
             [sind(MPCj.arrival_theta).*cosd(MPCj.arrival_phi),...
              sind(MPCj.arrival_theta).*sind(MPCj.arrival_phi),...
              cosd(MPCj.arrival_theta)]);
norm_MCD_AOA=sqrt(sum(MCD_AOA.^2,2));

MCD_AOD=0.5*([sind(MPCi.departure_theta).*cosd(MPCi.departure_phi),...
              sind(MPCi.departure_theta).*sind(MPCi.departure_phi),...
              cosd(MPCi.departure_theta)] - ...
             [sind(MPCj.departure_theta).*cosd(MPCj.departure_phi),...
              sind(MPCj.departure_theta).*sind(MPCj.departure_phi),...
              cosd(MPCj.departure_theta)]);
norm_MCD_AOD=sqrt(sum(MCD_AOD.^2,2));

delta_tau_max = max(MPCi.tau)-min(MPCi.tau);
tau_std = std(MPCi.tau);
norm_MCD_t = zeta*abs(MPCi.tau-MPCj.tau)/delta_tau_max*tau_std/delta_tau_max;

d = sqrt(norm_MCD_AOA.^2 + norm_MCD_AOD.^2 + norm_MCD_t.^2);

end

function meanAngle = meanAngles(angles)
% Author: Miead Tehrani-Moayyed
% Institute for the Wireless Internet of Things, 
% Northeastern University, Boston MA, 02115, USA
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 11-Sep-2022
%
% calculate the mean of a set of angles in degrees
%
% Input
%   angles    - a vector of angles in degrees
% Output
%   meanAngle - the mean of the angles 
%

angles = angles * pi/180;
angles = exp(1i*angles);
abar = mean(angles);
meanAngle = atan2(imag(abar),real(abar))*180/pi;
meanAngle(abs(abar)<1e-12) = nan;
end

function found = FindEmptyCluster(labels)
% Author: Miead Tehrani-Moayyed
% Institute for the Wireless Internet of Things, 
% Northeastern University, Boston MA, 02115, USA
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 11-Sep-2022
%
% This function checks if there are clusters with no member
% 
% Input
%   labels - a vector of labels
%
% Output
%   found - a flag to be set when empty clusters exist
%

found = 0;
for i = min(labels(:,1)):max(labels(:,1))
    nCk = sum(labels(:,1) == i);
    if  nCk == 0
        found = 1;
    end
end

end

function [tau_hat, h_hat] = FIR_coefficient_resampling(tau,h,nFIR,fs)
% arranges FIR coeefient from the input taps
%
% Input
%   tau - actual tap delays, tau (s)
%   h   - tap gain,
%   nFIR- Number of FIR filters
%   fs  - Sampling frequency
%
% Output
%   tau_hat - approximated taps delay of FIR filters
%   h_hat   - FIR Coefficients for approximated tap delay
%

if nargin >= 2
    if ~exist('nFIR','var')
        nFIR = 512; %number of FIR filters
    end
    if ~exist('fs','var')
        fs = 100e6; %100 M(Hz)
    end

end

Ts = 1/fs;
tau_hat = Ts:Ts:nFIR*Ts;

tapidx = round ((tau - 0) ./ Ts);

if sum(tapidx == 0) > 0     % dealing with the taps less than minimum excess delay, (e.g. 10ns for Colosseum)
    tapidx(tapidx == 0) = 1;       
    warning('taps exist with excess delay less than minimum excess delay')
end

if sum(tapidx > nFIR) > 0   % dealing with the taps greater than maximum excess delay (e.g. 5.12us for Colosseum)
    tapidx(tapidx > nFIR) = 512;
    warning('taps exist with excess delay greater than maximum excess delay')
end

if numel(tapidx) ~= numel(unique(tapidx))   % dealwing with the taps with similar excess delay
    warning('tap overlap accured in approximation')
end

h_hat = zeros(size(tau_hat));
for idx = 1:numel(h)
    h_hat(tapidx(idx)) = h_hat(tapidx(idx)) + h(idx);
end

tau_hat = tau_hat(unique(tapidx)).';
h_hat = h_hat(unique(tapidx)).';

end