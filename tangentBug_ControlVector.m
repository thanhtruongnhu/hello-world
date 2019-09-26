function u = tangentBug_ControlVector(id)
global Robot xp linesArr MOVE MAX MIN

New = zeros(1,2);

if Robot(id).mot2goal
    %% motion to xp(id,:)
    
    % update Robot(id).prevDist
    Robot(id).prevDist = Robot(id).currDist;
    
    % find distance&angle from Robot(id).x to xp(id,:)
    % invoke virtual pose sensors
    gl = getGoalPose(xp(id,:));
    rb = getRobotPose(Robot(id).x);
    
%     %  Visualization
%     figure(1);
%     plot(rb(1),rb(2),'g*');
%     hold on;
%     title('motion-to-xp(id,:)');
%     hold on;
    
    angRob2Goal = atan2(gl(2) - rb(2), gl(1) - rb(1));
    degAng = rad2deg(angRob2Goal);
    % read sensor
    reading = readRangeSensor(linesArr,rb); %linesArr: pairs of point for lines
    
    % check if there's obstacle on the direct line
    ang2Read = mod(round(degAng)+1,360); % Convert any angle into range [0;360]
    
    %% Check for obstacles on the line to xp(id,:)
    if  reading(ang2Read) == MAX % there's no obstacle
        Robot(id).goStr = 1;
    else % there's obstacle -> choose Oi to go straight to (instead of go straight to xp(id,:))
        Robot(id).goStr = 0;
        % we shall consider rob2Oi + Oi2goal as the distance function
        % instead of the direct distance from Robot(id).x to xp(id,:)
        if Robot(id).obsFol == 0 % This condition to prevent trigerring Boundary Following mode
            %                Cuz the Robot(id).prevDist in the last run Robot(id).step_for_fail_check (still in the state
            %                go straight to xp(id,:)) is the dist from Robot2Goal, this
            %                distance is definitely smaller than the new distance
            %                (which is calculated by rob2Oi + Oi2goal)
            Robot(id).prevDist = MAX;
            Robot(id).obsFol = 1;
        end
        
        % check if Robot(id).x is too close to obstacle
        if  reading(ang2Read) < MIN
            Robot(id).goTang = 1; % invoke tangent motion
        else
            Robot(id).goTang = 0;
        end
    end
    
    %% if Robot(id).x is following the path then motion-to-xp(id,:) is activated on
    % the straight line
    if Robot(id).goStr
        % follow straigth line to xp(id,:)
        New(1) = rb(1) + MOVE*cos(angRob2Goal);
        New(2) = rb(2) + MOVE*sin(angRob2Goal);
        Robot(id).currDist = findSqDistance(rb,gl);
        u = New - rb;
        if norm(u) ~=0
            u = (u/norm(u));
        end
        
        
        
    else % there's obstacle in front of Robot(id).x, find Oi's on the obstacle
        
        % find obstacles
        obs = findObstacles(reading,rb,gl);
        % a series of lines which demonstrate the discontinuity in [0;360]
        % This is a matrix with only 1 row
        
        % find the discontinuity points (pair of points)
        ind = find(obs >= round(degAng),1); % get the 1st element (angle) that greater than current heading
        if mod(ind,2) == 1 % odd location
            Oi = [obs(ind) , obs(ind+1)];
        else % even location
            Oi = [obs(ind-1) , obs(ind)];
        end
        
        if Oi(1) - Oi(2) == 0 % if the sensor sensed only a point   NEVER HAPPENS
            disp('WARNING:: Unexpected Behaviour !!! Check this !!!')
%             xNew = rb(1) + MOVE*cos(angRob2Goal);
%             yNew = rb(2) + MOVE*sin(angRob2Goal);
%             Robot(id).x(1) = xNew;
%             Robot(id).x(2) = yNew;
            Robot(id).currDist = findSqDistance(rb,gl);
            Robot(id).obsFol = 0;
        else % sensor has detected a segment intersecting the
            % line to the xp(id,:)
            
            % find the points of these disconts.
            %                 sensInd = mod(round(degAng)+1,360);
            P1 = [rb(1) + reading(Oi(1)+1)*cos(deg2rad(Oi(1))), ...
                rb(2) + reading(Oi(1)+1)*sin(deg2rad(Oi(1)))];
            
            P2 = [rb(1) + reading(Oi(2)+1)*cos(deg2rad(Oi(2))), ...
                rb(2) + reading(Oi(2)+1)*sin(deg2rad(Oi(2)))];
            
%             figure(1);
%             plot(P1(1),P1(2),'ko');
%             hold on;
%             plot(P2(1),P2(2),'ko');
%             hold on;
            
            % find the point that minimizes the rob2Oi + Oi2goal
            rob2O1 = findSqDistance(rb,P1);
            O12goal = findSqDistance(P1,gl);
            distO1 = rob2O1 + O12goal;
            
            rob2O2 = findSqDistance(rb,P2);
            O22goal = findSqDistance(P2,gl);
            distO2 = rob2O2 + O22goal;
            
            minOi = min(distO1,distO2);
            
            if Robot(id).prevPt > 0 % Robot(id).prevPt: selected Oi (This condition makes Robot(id).x keep going in the current direction)
                if Robot(id).prevPt == 1
                    minOi = distO1;
                else
                    minOi = distO2;
                end
            end
            
            % select the point that minimizes the distances
            if minOi == distO1 % select P1 (left side)
                angRob2Oi = atan2(P1(2)-rb(2), P1(1) - rb(1));
                Robot(id).prevPt = 1; % Robot(id).prevPt: selected Oi
                Robot(id).currDist = distO1;
            else % select P2 (right side)
                angRob2Oi = atan2(P2(2)-rb(2), P2(1) - rb(1));
                Robot(id).prevPt = 2; % Robot(id).prevPt: selected Oi
                Robot(id).currDist = distO2;
            end
            
            %% if distance started to increase initiate boundary following procedure
            if Robot(id).currDist > Robot(id).prevDist
                Robot(id).mot2goal = 0;
                Robot(id).ang2follow = angRob2Oi;
%                 Robot(id).dfollowed = minDist;
                Robot(id).firstBound = 1;
            else
                Robot(id).mot2goal = 1;
                
%                 % update min dist ever sensed
%                 if minOi == distO1
%                     minDist = O12goal;
%                 else
%                     minDist = O22goal;
%                 end
                
                
                if Robot(id).goTang % BEST PART OF TANGENT BUG (when Robot(id).x is too close to obstacle)
%                     degOi = rad2deg(findTangentAngle(angRob2Oi,Robot(id).prevPt));                                
                    New(1) = rb(1) + MOVE*cos(findTangentAngle(angRob2Oi, Robot(id).prevPt));
                    New(2) = rb(2) + MOVE*sin(findTangentAngle(angRob2Oi, Robot(id).prevPt));
                    u = New - rb;
                    if norm(u) ~=0
                        u = (u/norm(u));
                    end
                else % SHOULD REMOVE THIS PART, CUZ OUR SENSING RANGE IS VERY SMALL -> HENCE, OUR Robot(id).x COULD BUMP INTO OBSTACLE IF THIS CONDITION EXISTS
%                     degOi = rad2deg(angRob2Oi);                   
                    New(1) = rb(1) + MOVE*cos(angRob2Oi);
                    New(2) = rb(2) + MOVE*sin(angRob2Oi);
                    u = New - rb;
                    if norm(u) ~=0
                        u = (u/norm(u));
                    end
                end
            end
            
            
            
        end
        
    end
    
else
    %% boundary following
    % invoke virtual pose sensors
    gl = getGoalPose(xp(id,:));
    rb = getRobotPose(Robot(id).x);
    
    % angle to xp(id,:)
    angRob2Goal = atan2(gl(2)-rb(2), gl(1) - rb(1));
    degAng = rad2deg(angRob2Goal);
    
%     figure(1);
%     plot(rb(1),rb(2),'g*');
%     hold on;
%     title('boundary-following');
%     hold on;
    
    if Robot(id).firstBound
        Robot(id).boundStPos = rb; % Pose of Robot(id).x when boundary following mode started
    end
    
    % read sensor
    reading = readRangeSensor(linesArr,rb);
%     % check if there's obstacle on the direct line
%     ang2Read = mod(round(degAng)+1,360);
    
    obs = findObstacles(reading,rb,gl); % local angle (Robot(id).x)
    % find the discontinuity points
    ind = find(obs>=round(degAng),1);
    if mod(ind,2) == 1 % Odd number
        Oi = [obs(ind) , obs(ind+1)];
    else % even number
        Oi = [obs(ind-1) , obs(ind)];
    end
    % find the closest point to the Robot(id).x on the boundary
    [val,ind] =min(reading(Oi(1)+1 : Oi(2)+1));
    % reading: has numbering form 1:360 (MATLAB counting system)
    % Oi(i): has value ranging from 0:359 (Actually, findObstacles
    % returns value 0:360, but angle 0 is equay to 360)
    
    if val == 450 % saw the xp(id,:) ahead! or sonar max range
        % initiate motion to xp(id,:)
        Robot(id).mot2goal = 1;
        Robot(id).firstBound = 1;
        Robot(id).currDist = Robot(id).dreach;
        Robot(id).obsFol = 0;
        Robot(id).step_for_fail_check = 0;
    else % control Robot(id).x to move tangent with the closest point to it
        Prx = rb(1) + val*cos(deg2rad(Oi(1)+ind-1)); % -1 to return the right value of Oi(i)
        Pry = rb(2) + val*sin(deg2rad(Oi(1)+ind-1));
        Pr = [Prx Pry];
        
%         plot(Pr(1),Pr(2),'b*');
%         hold on;
        
        % distance between the xp(id,:) and the closest point on the boundary
        % to the Robot(id).x
        % angle between the closest pt and Robot(id).x
        angRob2Pr = atan2(Pr(2)-rb(2), Pr(1) - rb(1)); % global angle (global map)
        angRob2Pr = findTangentAngle(angRob2Pr,Robot(id).prevPt);
        
        [p,Robot(id).dreach] = find2ReachDist(reading,rb,gl); % Finding sensed point with minimum d(Oi,qgoal)
        
        if Robot(id).firstBound == 1 % Moving to  the Oi with nearest dist to xp(id,:)
            Robot(id).firstBound = 0;
            Robot(id).step_for_fail_check = Robot(id).step_for_fail_check + 1;
            Robot(id).minDFollowed = Robot(id).dfollowed;
            
            New(1) = rb(1) + MOVE*cos(Robot(id).ang2follow);
            New(2) = rb(2) + MOVE*sin(Robot(id).ang2follow);
            u = New - rb;
            if norm(u) ~=0
                u = (u/norm(u)); 
            end
            
        else % Boundary following mode
            
            if Robot(id).step_for_fail_check > 50 % Stop condition (Fail condition)
                if findSqDistance(Robot(id).boundStPos,rb) < 50 % after 50 run Robot(id).step_for_fail_check,
                    %        the distance btw the initial Robot(id).step_for_fail_check (when first starting boundary following
                    %        mode) and it's current position is less than 5 -> FAIL
                    Robot(id).failed = 1;
                end
            end
            
            [pt,Robot(id).dfollowed] = find2GoalDist(reading,Oi,rb,gl); % find point on obstacle boundary with minimum distance to xp(id,:)
            if Robot(id).dfollowed < Robot(id).minDFollowed
                Robot(id).minDFollowed = Robot(id).dfollowed;
            end
            % check for the end of boundary following
            if Robot(id).dreach < Robot(id).minDFollowed % switch back to boundary following
                Robot(id).mot2goal = 1;
                Robot(id).firstBound = 1;
                Robot(id).currDist = Robot(id).dreach;
                Robot(id).obsFol = 0;
                Robot(id).step_for_fail_check = 0;
            else                
                New(1) = rb(1) + MOVE*cos(angRob2Pr);
                New(2) = rb(2) + MOVE*sin(angRob2Pr);
                u = New - rb;
                if norm(u) ~=0
                    u = (u/norm(u));   
                end
            end
            
        end
        
    end
    
end


%% PLOTTING
% plot the Robot(id).x and covered path
%     figure(1)
% %     plot(xNew,yNew,'r*');
% %     hold on;
%     plot(Robot(id).x(1),Robot(id).x(2),'g.');
%     hold on;
%     if Robot(id).mot2goal
%         title('motion-to-xp(id,:)');
%     else
%         title('boundary-following');
%     end
%     hold on;
% plot the readings
%     reading = readRangeSensor(linesArr,Robot(id).x);
%     figure(2);
%     plot(reading);
%     axis([0,360,0,500]);
%     if Robot(id).failed
%         title('Robot(id).failed');
%         break;
%     end


%     if findSqDistance(Robot(id).x,xp(id,:)) < 5
%         title('SUCCESS');
%         break;
%     else
% %         clf(hnd);
%     end

end