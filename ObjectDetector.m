classdef ObjectDetector < matlab.System & matlab.system.mixin.CustomIcon & matlab.system.mixin.Propagates & matlab.system.mixin.SampleTime
    % OBJECTDETECTOR Object detector simulator
    %
    % Returns the angles, ranges, and labels of a simulated object detector
    %
    % For more information, see <a href="matlab:edit mrsDocObjectDetector">the documentation page</a>
    %
    % Copyright 2018-2019 The MathWorks, Inc.
    
    %% PROPERTIES
    % Public (user-visible) properties
    properties(Nontunable)
        mapName = ''; % Map
    end
    properties
        sensorOffset = [0,0];   % Object detector offset (x,y) [m] % coincidence with origin point of robot
        sensorAngle = 0;        % Object detector angle [rad] %Sensor orientation, in radians, relative to the vehicle x-axis 
        fieldOfView = pi/3;     % Sensor field of view [rad]
        maxRange = 5;           % Maximum range [m]
        maxDetections = 3;      % Maximum number of detections
    end
    properties(Nontunable)
       sampleTime = 0.1; % Sample time 
    end

    % Private properties
    properties(Access = private)
        map; % Occupancy grid
        hasMap;     % Binary flag for having a map
    end

    %% METHODS
    methods(Access = protected)
        
        % Setup method: Initializes all necessary graphics objects
        function setupImpl(obj)
            
            % Load the occupancy grid 
            obj.map = internal.createMapFromName(obj.mapName);
            obj.hasMap = ~isempty(obj.map);
            
        end

        % Step method: Outputs simulated object detections based on map,
        % robot pose, sensor field of view, max range, and object color
        function detections = stepImpl(obj,pose,objects)       
            
            % Initialize
            ranges = [];
            angles = [];
            labels = [];
            
            % Empty check
            if numel(objects) < 3
               detections = [];
               return; 
            end
            
            % Find the sensor pose and check if valid (check whether robot
            % is within display range)
            theta = pose(3);
            offsetVec = [cos(theta) -sin(theta); ...
                         sin(theta)  cos(theta)]*obj.sensorOffset'; %matrix with 2 rows, 2 columns
            sensorPose = pose' + [offsetVec', obj.sensorAngle]; % we don't have offset for sensor -> offsetVec matrix = 0, obj.sensorAngle = 0
            if ~obj.hasMap
                validPose = true;
            else
                validPose = sensorPose(1) >= obj.map.XWorldLimits(1) && ...
                            sensorPose(1) <= obj.map.XWorldLimits(2) && ... 
                            sensorPose(2) >= obj.map.YWorldLimits(1) && ... 
                            sensorPose(2) <= obj.map.YWorldLimits(2);
            end
                        
            % Return the RANGE and ANGLE for all objects
            % First, find the offsets
            offsets = objects(:,1:2) - sensorPose(1:2); % Find vector between object and sensor (in this case, sensor = robot)
            
            if ~isempty(offsets) && validPose
                % Extract RANGEs and ANGLEs
                ranges = sqrt(sum(offsets.^2,2)); % rows of matrix = DISTANCEs between object and sensor (in this case, sensor = robot)
                angles = wrapToPi(atan2(offsets(:,2),offsets(:,1))-sensorPose(3)); % rows of matrix = difference in angle 
%                 between object and sensor (in this case, sensor = robot)    
                
                % Filter by maximum range, field of vision**** -> can
                % either find impact angle and impact length on rod
                validIdx = (ranges <= obj.maxRange) & ... 
                           (abs(angles) <= obj.fieldOfView/2); % this will return a logical matrix (contains only 0 and 1)
                ranges = ranges(validIdx); % this will only get info of rows (with ID getting from validIdx) that contain 1 (True)
                angles = angles(validIdx); % this will only get info of rows (with ID getting from validIdx) that contain 1 (True)
                labels = objects(validIdx,3);               
             
           
                % Use occupancy grid, IF ANY, to account for obstacles
                % (maybe it don't run into this section)
                if ~isempty(obj.map)
                    % Loop backwards since we're removing values
                    for idx = numel(ranges):-1:1  
                        intPts = rayIntersection(obj.map,sensorPose, ... 
                                                 angles(idx),ranges(idx));
                        % Delete the reading if the point is occupied
                        if ~isnan(intPts)
                            ranges(idx) = [];
                            angles(idx) = [];
                            labels(idx) = [];
                        end
                    end
                end
                
            end
            
            if ~isempty(ranges)
                % Sort from nearest
                [ranges,sortedIdx] = sort(ranges);
                angles = angles(sortedIdx);
                labels = labels(sortedIdx);
                if numel(ranges) > obj.maxDetections
                    ranges = ranges(1:obj.maxDetections);
                    angles = angles(1:obj.maxDetections);
                    labels = labels(1:obj.maxDetections);
                end
                
                % Pack the final results into the output
                detections = [ranges, angles, labels];
            else
                detections = []; 
            end
                     
        end

        % More methods needed for the Simulink block to inherit its output
        % sizes from the scan angle parameter provided.
        function sz = getOutputSizeImpl(obj)
            sz = [obj.maxDetections,3];
        end
        
        function fx = isOutputFixedSizeImpl(~)
           fx = false;
        end
        
        function dt = getOutputDataTypeImpl(obj)
            dt = propagatedInputDataType(obj,1);
        end

        function cp = isOutputComplexImpl(~)
            cp = false;
        end
        
        % Define icon for System block
        function icon = getIconImpl(~)
            icon = {'Object','Detector'};
        end
        
        % Define sample time
        % Must be discrete since the output is a variable-sized signal
        function sts = getSampleTimeImpl(obj)
            if obj.sampleTime > 0
               sts = createSampleTime(obj,'Type','Discrete','SampleTime',obj.sampleTime);
            else
               sts = createSampleTime(obj,'Type','Inherited'); 
            end
        end
        
        % Save and load object methods to enable stepping
        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            s.map = obj.map;
        end
        function loadObjectImpl(obj,s,wasInUse)
            obj.map = s.map;
            loadObjectImpl@matlab.System(obj,s,wasInUse);
        end
        
    end
    
    methods (Static, Access = protected)
        % Do not show "Simulate using" option
        function flag = showSimulateUsingImpl
            flag = false;
        end
        % Always run in interpreted mode
        function simMode = getSimulateUsingImpl
            simMode = 'Interpreted execution';
        end
    end
    
end
