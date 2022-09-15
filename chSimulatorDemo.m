% Ray-tracing Based Mobile Wireless Channel Simulator Demo
% Mobile Scenario Example - NU campus LTE small cell and WIFI coexistence
%
% Author: Miead Tehrani-Moayyed
% Institute for the Wireless Internet of Things, 
% Northeastern University, Boston MA, 02115, USA
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 11-Sep-2022

% Parameters Configuration
parameters.origin = [42.34025, -71.08848];              %Wireless environment origin
parameters.osmfile = "NUcampus.osm";                    %Wireless environment 3D Model OSM file
parameters.fc = 5.8e9;                                  %LTE ISM band Center freq. (Hz)
parameters.nodes.names = {'LTE' 'WIFI' 'UE#1' 'UE#2' 'WIFIUE'};  %Nodes names
parameters.nodes.powers = {30, 22, 20, 20, 15};         %Transmit power (dBm)
parameters.nodes.coordinates = {[42.3405645965253,-71.0892925339939],...        %stationary nodes coordinates
                                [42.3398095978303,-71.0884277156665],...
                                [42.341482,-71.086650; 42.341102,-71.087238; 42.339903,-71.090416],...
                                [42.339757,-71.090343; 42.341045,-71.087006; 42.341412,-71.086518]};%Control points for route
parameters.nodes.heights = {17, 17, [1.5, 1.5, 1.5], [1.5, 1.5, 1.5], 1.5};                         %Height of nodes
parameters.nodes.antenna = {'isotropic' 'isotropic' 'isotropic' 'isotropic' 'isotropic'};           %Nodes antenna type
parameters.nodes.antAngle = {0 0 0 0 0};                %Nodes antenna azimuth angles 
parameters.nodes.velocity = {0 0 11.18 5.6 nan};        %Nodes velocity [m/s]
parameters.nodes.mobilityType = {'stationary' 'stationary' 'route' 'route' 'custom'};
parameters.Ts = 1000e-3;                                %Mobile channel sampling interval [s]
parameters.scenarioDuration = 70;                       %Scenario duration [s]
parameters.rt.nReflections = 3;                         %Number of reflection in raytracing
parameters.rt.BuildingsMaterial = "custom";             %Type of buildings material
parameters.rt.BuildingsMaterialPermittivity = 5.31;     %Building material permitivity for custom material
parameters.rt.BuildingsMaterialConductivity = 0.1353;   %Building material conductivity for custom material
parameters.rt.TerrainMaterial = "custom";               %Type of terrain material
parameters.rt.TerrainMaterialPermittivity = 5.31;       %Terrain material permitivity for custom material
parameters.rt.TerrainMaterialConductivity = 0.1353;     %Terrain material conductivity for custom material
parameters.rt.polarization = "H";                       %Polarization
parameters.visualization.flag = 1;                      %Visualization request
parameters.visualization.TxId = 1;                      %TX Index for visualization
parameters.visualization.RxId = 4;                      %RX Index for visualization

% Define the custom mobility pattern - RWP in this case
RWP = load('RWPmobility.mat');
parameters.nodes.coordinates {5} = [RWP.snapshots.lat RWP.snapshots.lon]; 

% Channel simulation process
chMatrix = channelSimulator(parameters);

% Plot time evolution of paths
toa_scale_factor = 1e6;

figure
TxID = parameters.visualization.TxId;
RxID = parameters.visualization.RxId;
for snapshotIdx = 1:size(chMatrix,3)
    scatter( chMatrix{TxID,RxID,snapshotIdx}.tau * toa_scale_factor, ...
             ones(size(chMatrix{TxID,RxID,snapshotIdx}.tau)) * parameters.Ts * snapshotIdx, [], ...
             mag2db(abs(chMatrix{TxID,RxID,snapshotIdx}.h)) , 'filled');
    hold on
end
colorbar
grid on
xlabel('TOA [$\mu S$]','Interpreter','latex');
ylabel('Time [s]')

% Plot Path Loss and Generate path loss matrix
TxID = [1, 2];
RxID = 4;

plMatrix = generatePLmatrix(chMatrix);

% plot Path Loss using plMatrix
figure
legends = strings(numel(TxID)*numel(RxID),1);
for TxIdx = 1 : numel(TxID)

    for RxIdx = 1 : numel(RxID)
        Tx = (TxID(TxIdx));
        Rx = (RxID(RxIdx));

        PL = plMatrix(Tx,Rx,:);

        plot(Time/1000,PL(:));
        hold on
        legends((TxIdx-1)*numel(RxID)+RxIdx) = sprintf("%s to %s",parameters.nodes.names{Tx},parameters.nodes.names{Rx});

    end

end
legend(legends);
ylabel("Path Loss [dB]");
xlabel("Time [s]");
grid on;
ylim([50,150])

% Channel taps approximation
% Approximate the channels to the four taps for emulation on Colosseum
[chMatrix_hat, plMatrix_hat] = chApproximationV1(chMatrix,'elevation','NLOS');

% Plot time evolution of approximated taps
TxID = 1;
RxID = 4;
toa_scale_factor=1e6;

figure(6)
tau =[];
for snapshotIdx = 1:size(chMatrix,3)
    scatter( cell2mat(chMatrix_hat{TxID,RxID,snapshotIdx}(:,1)) * toa_scale_factor, ...
             ones(size(chMatrix_hat{TxID,RxID,snapshotIdx}(:,1))) * parameters.Ts * snapshotIdx, [], ...
             mag2db(abs(cell2mat(chMatrix_hat{TxID,RxID,snapshotIdx}(:,2)))) , 'filled');
    hold on
end
colorbar
grid on
xlabel('TOA [$\mu S$]','Interpreter','latex');
ylabel('Time [s]')