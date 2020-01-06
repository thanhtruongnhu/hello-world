function ControlVector_3(r_idx)
% Level 1: Behavioural Control
% Calculate Ouput Velocity Vector [Vx,Vy] for Robot(r_idx)

global Tc epsilon r_a
global Robot index

Vs = zeros(1,2); % only in Fuzzy Controller
Vc = zeros(1,2);
Vc_fuzzy = zeros(1,2);
Robot(r_idx).Vi_max = (epsilon/2)/Tc;
gamma = 1;


if ~isempty(Robot(r_idx).nrst_pnt) %Detect something in obstacle avoidance zone
    %% 1. Fuzzy Controller (Wall Following Behavior - Obstacle avoidance zone)
    index = r_idx; % [Visualization] Input data for Visualize Fuzzy Data button    
    
    r_ij = Robot(r_idx).nrst_pnt - Robot(r_idx).x;
    r = norm(r_ij);
    r_unit = r_ij/r;
    Vs_comp_r2o = -r_unit;
    Vc_comp_r2o = r_unit;
    Vs = Vs + Vs_comp_r2o
    Vc_fuzzy = Vc_fuzzy + Vc_comp_r2o;    
    
    % Va = allignVel_3(r_idx, Robot(r_idx).nrst_pnt);    % Approach 1: Create a vector perpendicular to vector [robot;obstacle]
    % Va = allignVel_4(r_idx) % Approach 2: V_a keep robot going in the direction of its current heading
    Va = allignVel_5(r_idx, Robot(r_idx).nrst_pnt); % Approach 3: Create a vector perpendicular to vector [robot;obstacle] but steer the robot to the right side only !
    Robot(r_idx).v_a = Va;
    
    [alpha,beta] = Fuzzy_controller(r_idx); 
    Robot(r_idx).v = alpha*Vc_fuzzy + beta*Vs + gamma*Va;
else %Detect nothing in obstacle avoidance zone
    %% 2. Pure Behavior Controller (Free zone & Critical zone)
    if ~isempty(Robot(r_idx).neighbor)
        for k = 1 : size(Robot(r_idx).neighbor, 1)
%             size(Robot(r_idx).neighbor, 1)
%             Robot(r_idx).neighbor(k,2:3)

            if r > r_a  % obstacle avoidance zone radius (r_a) = tactile sensing radius (MAXREADING)
                r_ij = Robot(r_idx).neighbor(k,2:3) - Robot(r_idx).x;
                %         [Vs_comp_r2r,Vc_comp] = SeparateAndCoherent(r_ij);
                r = norm(r_ij);
                r_unit = r_ij/r;
                Vc_comp = r_unit;
                Vc = Vc + Vc_comp;  % Vector combination (Coherent Velocity)
            end
        end
        
    else
        %     Vs = [0 0];
        Vc = [0 0];
    end
    Robot(r_idx).v = 1*Vc; % temporary weight = 1
    % Va = 0.5*allignVel_2(r_idx, xdk);       % Allignment Velocity
    % Robot(r_idx).v = gamma*(sigma*Vc + Va) + rho*Vs; % CONTROL VECTOR
end
V_observe = Robot(r_idx).v  %observer only


%% Max velocity constraint: v_i_max = epsilon/2/delta(t)
if norm(Robot(r_idx).v) > Robot(r_idx).Vi_max
    Robot(r_idx).v = Robot(r_idx).Vi_max*Robot(r_idx).v/norm(Robot(r_idx).v) ;
end


end


