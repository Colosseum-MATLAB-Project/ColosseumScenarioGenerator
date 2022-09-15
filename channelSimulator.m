function [chMatrix, coordinates, Time]= channelSimulator(parameters)
% Channel simulator function
%
% Author: Miead Tehrani-Moayyed
% Institute for the Wireless Internet of Things, 
% Northeastern University, Boston MA, 02115, USA
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 11-Sep-2022
%
% Simulate mobile wireless channels using ray-tracing propagation model
% For visualization, set function visualization input parameter
%
% Input: paramaters as an struct with the following format
% parameters.origin                             Wireless environment origin
% parameters.osmfile                            Wireless environment 3D Model OSM file
% parameters.fc                                 Carrier frequency (Hz)
% parameters.nodes.names                        Nodes names cell array
% parameters.nodes.powers                       Transmit power (dBm)
% parameters.nodes.coordinates                  Nodes coordinates/control points for stationary/mobile 
% parameters.nodes.heights                      Height of nodes
% parameters.nodes.antenna                      Nodes antenna type
% parameters.nodes.antAngle                     Nodes antenna azimuth angles 
% parameters.nodes.velocity                     Nodes velocity [m/s] scalar values cell 
% parameters.nodes.mobilityType                 Nodes mobility type
% parameters.Ts                                 Mobile channel sampling interval [s]
% parameters.scenarioDuration                   Scenario duration [s]
% parameters.rt.nReflections                    Number of reflections in raytracing
% parameters.rt.BuildingsMaterial               Type of buildings material
% parameters.rt.BuildingsMaterialPermittivity   Building material permitivity for custom material
% parameters.rt.BuildingsMaterialConductivity   Building material conductivity for custom material
% parameters.rt.TerrainMaterial                 Type of terrain material
% parameters.rt.TerrainMaterialPermittivity     Terrain material permitivity for custom material
% parameters.rt.TerrainMaterialConductivity     Terrain material conductivity for custom material
% parameters.rt.polarization                    Polarization {"H" or "V"}
% parameters.visualization.flag                 Visualization flag {1 or 0}
% parameters.visualization.TxId                 TX Index for visualization
% parameters.visualization.RxId                 RX Index for visualization
%
% Output:
% chMatrix                                      Channel matrix table
%
% Example: see chSimulatorDemo.mlx or chSimulatorDemo.m

% Obtaining samples for mobile nodes
positions = cell(size(parameters.nodes.names));
for nodeIdx = 1 : numel(parameters.nodes.names)
    % Obtaining the coordinates for mobile nodes defined by route
    if strcmp(parameters.nodes.mobilityType{nodeIdx},'route') == 1

        spacing = parameters.nodes.velocity{nodeIdx} * parameters.Ts;
        [points, ~] = createRoute('CoordinateSystem','geographic',...
                     'Latitudes', parameters.nodes.coordinates{nodeIdx}(:,1),...
                     'Longitudes',parameters.nodes.coordinates{nodeIdx}(:,2),...
                     'Heights', parameters.nodes.heights{nodeIdx},...
                     'Spacing',spacing);
        positions{nodeIdx} = [points.lats' points.lons'];

    % Coordinates for custom mobile nodes
    elseif strcmp(parameters.nodes.mobilityType{nodeIdx},'custom') == 1
        positions{nodeIdx} = parameters.nodes.coordinates{nodeIdx};
    
    % Coordinates for stationary nodes    
    elseif strcmp(parameters.nodes.mobilityType{nodeIdx},'stationary') == 1
        positions{nodeIdx} = parameters.nodes.coordinates{nodeIdx};

    end
end

% Generating snapshots for mobility simulation
% Obtaining the number of snapshots for the duration of scenario
nSnapshots = parameters.scenarioDuration / parameters.Ts;

% Constructing nodes coordinates for snapshots
coordinates = nan(numel(parameters.nodes.names),3,nSnapshots);
for snapshotIdx = 1:nSnapshots
    for nodeIdx = 1:numel(parameters.nodes.names)
        if parameters.nodes.velocity{nodeIdx} == 0
            coordinates(nodeIdx,:,snapshotIdx) = [positions{nodeIdx} parameters.nodes.heights{nodeIdx}];
        else
            coordinates(nodeIdx,:,snapshotIdx) = [positions{nodeIdx}(min(snapshotIdx,size(positions{nodeIdx},1)),:) parameters.nodes.heights{nodeIdx}(1)];
        end
    end

end

Time = (0:nSnapshots-1) * parameters.Ts*1e3;

fprintf('%d snapshots generated for simulation\n',nSnapshots);

% Defining wireless environment and nodes 
% Loading wireless environment
viewer = siteviewer("Buildings",parameters.osmfile,"Basemap","satellite");


% Defining Transmitting nodes
txs = txsite('Name', parameters.nodes.names,...
       'TransmitterPower',db2pow(cell2mat(parameters.nodes.powers))/1000, ...
       'Latitude',coordinates(:,1,1),...
       'Longitude',coordinates(:,2,1), ...
       'Antenna',parameters.nodes.antenna,...
       'AntennaAngle',cell2mat(parameters.nodes.antAngle),...
       'AntennaHeight',coordinates(:,3,1)',...
       'TransmitterFrequency',parameters.fc);
   
% Defining Receiving nodes
rxs = rxsite('Name', parameters.nodes.names,...
       'Latitude',coordinates(:,1,1),...
       'Longitude',coordinates(:,2,1), ...
       'Antenna',parameters.nodes.antenna,...
       'AntennaAngle',cell2mat(parameters.nodes.antAngle),...
       "AntennaHeight",coordinates(:,3,1)');

% Mobile Scenario Visualization
snapshot = 1;       %Determine the snapshot to be visualized, defualt: 1

visualization = parameters.visualization.flag;
if visualization
    clearMap(viewer);
    for nodeIdx = 1:numel(txs)
        pause(5)

        txs(nodeIdx).Latitude = coordinates(nodeIdx,1,snapshot);
        txs(nodeIdx).Longitude = coordinates(nodeIdx,2,snapshot);
        txs(nodeIdx).AntennaHeight = coordinates(nodeIdx,3,snapshot);

        show(txs(nodeIdx))

        %Show antenna pattern other than isotropic
        if ~strcmp(txs(nodeIdx).Antenna , 'isotropic')
            pause(5)
            pattern(txs(nodeIdx))
        end

        % show samples of the mobile nodes
        if parameters.nodes.velocity{nodeIdx} ~= 0
            points = txsite('Name', strcat(repmat(strcat(parameters.nodes.names(nodeIdx), '-'),1,nSnapshots),arrayfun(@num2str, 1:nSnapshots, 'UniformOutput', 0) ), ...
                'Latitude',reshape(coordinates(nodeIdx,1,:),[1,nSnapshots]),...
                'Longitude',reshape(coordinates(nodeIdx,2,:),[1,nSnapshots]), ...
                "AntennaHeight",reshape(coordinates(nodeIdx,3,:),[1,nSnapshots]));

            show(points,...
                'Icon',"Gdot.png",...
                "IconSize",[5 5]), pause(10)
        end

    end
end

% Mobile node snapshots animation
% Show mobile node snapshots animation
ffRatio = 32;               % Fast Forward ratio

if visualization

    for snapshotIdx = 1:nSnapshots

        clearMap(viewer);

        points = txsite('Name', parameters.nodes.names, ...
            'Latitude',coordinates(:,1,snapshotIdx),...
            'Longitude',coordinates(:,2,snapshotIdx), ...
            "AntennaHeight",coordinates(:,3,snapshotIdx)');
        show(points)

        pause(parameters.Ts/ffRatio);

    end

end

% Raytracing simulation
% Define propagation prediction model
rtpm = propagationModel("raytracing", ...
    "Method","sbr", ...
    "AngularSeparation","low", ...
    "MaxNumReflections",parameters.rt.nReflections, ...
    "BuildingsMaterial",parameters.rt.BuildingsMaterial, ...
    "BuildingsMaterialPermittivity", parameters.rt.BuildingsMaterialPermittivity, ...
    "BuildingsMaterialConductivity", parameters.rt.BuildingsMaterialConductivity, ...
    "TerrainMaterial",parameters.rt.TerrainMaterial, ...
    "TerrainMaterialPermittivity", parameters.rt.TerrainMaterialPermittivity, ...
    "TerrainMaterialConductivity",parameters.rt.TerrainMaterialConductivity);

disp('Ray-tracing snapshots...');
chData = cell(numel(txs),numel(rxs),nSnapshots);
tic
for snapshotIdx =  1:nSnapshots

    %update location of nodes
    for nodeIdx = 1:numel(txs)
        txs(nodeIdx).Latitude = coordinates(nodeIdx,1,snapshotIdx);
        txs(nodeIdx).Longitude = coordinates(nodeIdx,2,snapshotIdx);
        txs(nodeIdx).AntennaHeight = coordinates(nodeIdx,3,snapshotIdx);

        rxs(nodeIdx).Latitude = coordinates(nodeIdx,1,snapshotIdx);
        rxs(nodeIdx).Longitude = coordinates(nodeIdx,2,snapshotIdx);
        rxs(nodeIdx).AntennaHeight = coordinates(nodeIdx,3,snapshotIdx);
    end

    chData(:,:,snapshotIdx) = raytrace(txs,rxs,rtpm);

    fprintf('%d,',snapshotIdx);
end
fprintf('\n');
toc

fprintf('\n');

%Apply polarization for ray-traced channels
Polarization = parameters.rt.polarization;

for snapshotIdx = 1:size(chData,3)
    for txidx = 1 : size(chData,1)
        for rxidx = 1 : size(chData,2)
            if txidx==rxidx, continue, end
            for rayidx = 1 : numel(chData{txidx,rxidx,snapshotIdx})
                [chData{txidx,rxidx,snapshotIdx}(rayidx).PathLoss,chData{txidx,rxidx,snapshotIdx}(rayidx).PhaseShift] = raypl(chData{txidx,rxidx,snapshotIdx}(rayidx),...
                    "TransmitterPolarization",Polarization,...
                    "ReceiverPolarization",Polarization);
            end
        end
    end
end

% Generate channel matrix
chMatrix = cell(size(chData));
for snapshotIdx = 1:size(chData,3)
    for txIdx = 1:size(chData,1)
        for rxIdx = 1:size(chData,2)

            nRays = length(chData{txIdx,rxIdx,snapshotIdx});
            rays = struct (...
                'tau',cell(nRays,1),...
                'h',cell(nRays,1),...
                'arrival_theta',cell(nRays,1),...
                'arrival_phi',cell(nRays,1),...
                'departure_theta',cell(nRays,1),...
                'departure_phi',cell(nRays,1));

            for rayIdx = 1:length(chData{txIdx,rxIdx,snapshotIdx})
                rays(rayIdx).h = db2mag(-chData{txIdx,rxIdx,snapshotIdx}(1,rayIdx).PathLoss) .* exp(chData{txIdx,rxIdx,snapshotIdx}(1,rayIdx).PhaseShift*1i);
                rays(rayIdx).tau = chData{txIdx,rxIdx,snapshotIdx}(1,rayIdx).PropagationDelay;
                rays(rayIdx).departure_phi = chData{txIdx,rxIdx,snapshotIdx}(1,rayIdx).AngleOfDeparture(1);
                rays(rayIdx,1).departure_theta = chData{txIdx,rxIdx,snapshotIdx}(1,rayIdx).AngleOfDeparture(2);
                rays(rayIdx,1).arrival_phi = chData{txIdx,rxIdx,snapshotIdx}(1,rayIdx).AngleOfArrival(1);
                rays(rayIdx,1).arrival_theta = chData{txIdx,rxIdx,snapshotIdx}(1,rayIdx).AngleOfArrival(2);
            end

            chMatrix{txIdx,rxIdx,snapshotIdx} = struct2table(rays,'AsArray',true);

        end

    end
end

% Visualize propagation path and time-varing CIR
TxID = parameters.visualization.TxId;
RxID = parameters.visualization.RxId;
toa_scale_factor=1e6;
colors = ['b'; 'r';'g';'m'; 'c'; 'k'; 'b'; 'y'];
nPaths = 20;

for snapshotIdx = 1 : size(chMatrix,3)
    figure(5)
    %set(gcf,'Visible','on')    % for external figure for recording
    legends = strings(1,1);
    legendsIdx = 1;
    for TxIdx = 1 : numel(TxID)
        for RxIdx = 1 : numel(RxID)
            Tx = TxID(TxIdx);
            Rx = RxID(RxIdx);

            if size(chMatrix{Tx,Rx,snapshotIdx},1) > 0
                
                    stem((chMatrix{Tx,Rx,snapshotIdx}.tau)*toa_scale_factor, ...
                        mag2db(abs(chMatrix{Tx,Rx,snapshotIdx}.h)),'BaseValue',(-150), ...
                        'Color',colors((TxIdx-1)*numel(RxID)+RxIdx),'LineWidth',1.25);

                    legends(legendsIdx) = sprintf("%s to %s",parameters.nodes.names{Tx},parameters.nodes.names{Rx});
                    legend(legends);
                    legendsIdx = legendsIdx + 1;
               
            else
                stem([],[]);
            end

            hold on
        end
    end

    hold off
    ylim([-150,-60])
    xlim([0,3.5])
    ylabel('Path Gain [dB]');
    xlabel('TOA [$$\mu S$$]','Interpreter','latex');
    grid on
    text(xlim*[7/8;1/8],ylim*[1/8;7/8]+1,sprintf('Time: %.2f [S]',Time(snapshotIdx)*1e-3));

    if visualization
    if size(chData{1,snapshotIdx},2)>0
        clearMap(viewer)

        txs(TxID(1)).Latitude = coordinates(TxID(1),1,snapshotIdx);
        txs(TxID(1)).Longitude = coordinates(TxID(1),2,snapshotIdx);
        txs(TxID(1)).AntennaHeight = coordinates(TxID(1),3,snapshotIdx);

        rxs(RxID(1)).Latitude = coordinates(RxID(1),1,snapshotIdx);
        rxs(RxID(1)).Longitude = coordinates(RxID(1),2,snapshotIdx);
        rxs(RxID(1)).AntennaHeight = coordinates(RxID(1),3,snapshotIdx);

        show(txs(Tx),'ShowAntennaHeight', true)
        show(rxs(Rx),'ShowAntennaHeight', true)

        rays = chData{TxID(1),RxID(1),snapshotIdx};
        plot(rays(:,1:min(size(rays,2),nPaths)),'Type', 'pathloss');
    end
    end
    ffRatio = 16;               % Fast Forward ratio
    pause(parameters.Ts/ffRatio);   
    
    

end

end