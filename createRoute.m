function [position, names] = createRoute (varargin) 
% Project Name: Creating Route for outdoor and indoor scenario  
%               
% File Name: createRoute.m
%
% Author: Miead Tehrani-Moayyed
% Work address: Wireless Networks and System Lab  
% Northeastern University, 360 Huntington Ave. Boston, MA 02115
% email: tehranimoayyed.m@northeastern.edu
% Last revision: 22-Apr-2022
%
% This function creates a route using control points as the
% points to form the trajectory.
%
% Properties:
%   
%   CoordinateSystem: 'geographic' for outdoor or 'cartesian' for indoor
%   Latitudes       - latitudes array for the control points
%   Longitudes      - longitudes array the control points
%   Positions:      - cartesian coordinates for indoor control points defines as
%                     a 3 by n matrix with the format [x1,x2..;y1,y2..;z1,z2..]
%   Spacing:        - route element spacing
%
% Outputs
%   position         
%                   - For geographic coordinates
%                       position.lats: route points latitudes
%                       position.lons: route points longitudes
%                   - For cartesian coordinates includes [x;y;z]
%   names           - route points name including the location tag number
%
%   Examples
%   Example 1: Creating a route in an outdoor scenario
%        [positions, names] = createRoute('CoordinateSystem','geographic',...
%                             'Latitudes', [42.340257, 42.340487, 42.339876, 42.340100],...
%                             'Longitudes',[-71.088774, -71.088241, -71.088474, -71.087960 ],...
%                             'Heights', [2, 2, 2, 2],...
%                             'Spacing',7);
%
%
%   Example 2: Creating a route in an indoor scenario
%         [positions, names] = createRoute ('CoordinateSystem','cartesian',...
%                             'Positions', [-1,1,-1,1;1,1,-1,-1;.85,.85,0.85,.85],...
%                               'Spacing', 0.25);
%
%
% ------------- BEGIN CODE --------------

defaultRouteLatitudes = [42.340257, 42.340487, 42.339876, 42.340100];
defaultRouteLongitudes = [-71.088774, -71.088241, -71.088474, -71.087960 ];
defaultRouteHeights = [2, 2, 2, 2];
defaultIndoorPosition = [-1,1,-1,1;1,1,-1,-1;.85,.85,0.85,.85];

p = inputParser;
validGeoCtrlPts = @(x) isnumeric(x) && isvector(x) && (numel(x) > 2);
addOptional(p,'CoordinateSystem','geographic', @(x) ischar(x) && (strcmp(x,'geographic') || strcmp(x, 'cartesian')))
addOptional(p,'Latitudes', defaultRouteLatitudes, validGeoCtrlPts )  % Quadrant in NU campus as default
addOptional(p,'Longitudes', defaultRouteLongitudes, validGeoCtrlPts)   
addOptional(p,'Heights', defaultRouteHeights, validGeoCtrlPts)   
addOptional(p,'Positions', defaultIndoorPosition, @(x) isnumeric(x) && ismatrix(x) && sum(size(x) == [3,4] )==2 ), 
addOptional(p,'Spacing',1, @(x) isscalar(x) && (x>0));

parse(p,varargin{:});

route.controlPoints.coordinateSystem = p.Results.CoordinateSystem;
route.controlPoints.lats = p.Results.Latitudes;
route.controlPoints.lons = p.Results.Longitudes;
route.controlPoints.heights = p.Results.Heights;
route.controlPoints.positions = p.Results.Positions;
route.spacing = p.Results.Spacing;

if strcmpi(route.controlPoints.coordinateSystem,'geographic')
    % Define route control points as the sites
    ctrlPts = rxsite("Latitude",route.controlPoints.lats, ...
        "Longitude",route.controlPoints.lons, ...
        "AntennaHeight",route.controlPoints.heights);

    % Calculate parameters for obtaining route point location
    azs = nan(1,numel(ctrlPts)-1);
    nSensor = nan(1,numel(ctrlPts)-1);
    for ptsIdx = 1 : numel(ctrlPts)-1
        azs(ptsIdx) = angle(ctrlPts(ptsIdx),ctrlPts(ptsIdx+1));
        dist = distance(ctrlPts(ptsIdx),ctrlPts(ptsIdx+1));
        nSensor(ptsIdx) = floor(dist/route.spacing);
    end

    % Calculate route element locations
    lats = nan(1,sum(nSensor));
    lons = nan(1,sum(nSensor));
    names = strings(1,sum(nSensor));
    idx = 1;
    for ptsIdx = 1 : numel(ctrlPts)-1
        for sensorIdx = 1 : nSensor(ptsIdx)
            d = (sensorIdx-1) * route.spacing;
            az = azs(ptsIdx);

            [lats(idx),lons(idx)] = location(ctrlPts(ptsIdx),d,az);
            names(idx) = sprintf(" P%d",idx);

            idx = idx +1;
        end

    end
    
    position.lats = lats;
    position.lons =lons;

elseif strcmpi(route.controlPoints.coordinateSystem,'cartesian')
    % Define route control points as the indoor sites
    ctrlPts = rxsite("cartesian","AntennaPosition",route.controlPoints.positions);

    % Calculate parameters for obtaining route point location
    azs = nan(1,numel(ctrlPts)-1);
    nSensor = nan(1,numel(ctrlPts)-1);
    for ptsIdx = 1 : numel(ctrlPts)-1
        azs(ptsIdx) = angle(ctrlPts(ptsIdx),ctrlPts(ptsIdx+1));
        dist = distance(ctrlPts(ptsIdx),ctrlPts(ptsIdx+1));
        nSensor(ptsIdx) = floor(dist/route.spacing);
    end

    % Calculate route element locations
    Xs = nan(1,sum(nSensor));
    Ys = nan(1,sum(nSensor));
    Zs = nan(1,sum(nSensor));
    names = strings(1,sum(nSensor));
    idx = 1;
    for ptsIdx = 1 : numel(ctrlPts)-1
        for sensorIdx = 1 : nSensor(ptsIdx)
            d = (sensorIdx-1) * route.spacing;
            az = azs(ptsIdx);

            Xs(idx) = ctrlPts(ptsIdx).AntennaPosition(1) + d * cosd(az);
            Ys(idx) = ctrlPts(ptsIdx).AntennaPosition(2) + d * sind(az);
            Zs(idx) = ctrlPts(ptsIdx).AntennaPosition(3);
            names(idx) = sprintf(" P%d",idx);

            idx = idx +1;
        end

    end

    position = [Xs;Ys;Zs];
end

end