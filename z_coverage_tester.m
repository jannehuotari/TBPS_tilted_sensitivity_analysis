%% z-coverage

%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


layer = 1; % Layer to analyze
modules = 'outer'; % "inner" or "outer" modules to analyze
changes = 'no'; % Include random variation to module r, z and tilt angle (yes/no)
tolerances = 'no'; % Calculate min / max tolerance values for module displacement (yes/no)
runs = 1; % Simulation runs
accuracy = 0.01; % Step size of tolerance calculation
sensor_thickness = 0.200; % Thickness of the sensors
tilt_var_min = 0; % Variance min of module tilt angle
tilt_var_max = 0; % Variance max of module tilt angle
r_var_min = 0; % Variance min of module position in r
r_var_max = 0; % Variance max of module position in r
z_var_min = 0; % Variance min of module position in z
z_var_max = 0; % Variance max of module position in z


%%%%%%%%%%%%%%%%%%%%%%%%% Loop Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize some arrays
z_cov = zeros(1,12);

r_results = zeros(runs,12);
tilt_results = zeros(runs,12);
z_results = zeros(runs,12);
z_cov_results = zeros(runs,12);

r_s1 = zeros(1,12);
r_s2 = zeros(1,12);
tilt_m = zeros(1,12);
z_s1 = zeros(1,12);
z_s2 = zeros(1,12);

z_tolerance = zeros(2,13);

% Nominal module positions:
if layer == 1
    if strcmp('inner', modules) == 1
        r_s1_i = r_s1_layer1_inner;
        r_s2_i = r_s2_layer1_inner;
        tilt_m_i = tilt_layer1;
        z_s1_i = z_s1_layer1_inner;
        z_s2_i = z_s2_layer1_inner;
    elseif strcmp('outer', modules) == 1
        r_s1_i = r_s1_layer1_outer;
        r_s2_i = r_s2_layer1_outer;
        tilt_m_i = tilt_layer1;
        z_s1_i = z_s1_layer1_outer;
        z_s2_i = z_s2_layer1_outer;
    end
elseif layer == 2
    if strcmp('inner', modules) == 1
        r_s1_i = r_s1_layer2_inner;
        r_s2_i = r_s2_layer2_inner;
        tilt_m_i = tilt_layer2;
        z_s1_i = z_s1_layer2_inner;
        z_s2_i = z_s2_layer2_inner;
    elseif strcmp('outer', modules) == 1
        r_s1_i = r_s1_layer2_outer;
        r_s2_i = r_s2_layer2_outer;
        tilt_m_i = tilt_layer2;
        z_s1_i = z_s1_layer2_outer;
        z_s2_i = z_s2_layer2_outer;
    end
elseif layer == 3
    if strcmp('inner', modules) == 1
        r_s1_i = r_s1_layer3_inner;
        r_s2_i = r_s2_layer3_inner;
        tilt_m_i = tilt_layer3;
        z_s1_i = z_s1_layer3_inner;
        z_s2_i = z_s2_layer3_inner;
    elseif strcmp('outer', modules) == 1
        r_s1_i = r_s1_layer3_outer;
        r_s2_i = r_s2_layer3_outer;
        tilt_m_i = tilt_layer3;
        z_s1_i = z_s1_layer3_outer;
        z_s2_i = z_s2_layer3_outer;
    end
end

r_s1 = r_s1_i;
r_s2 = r_s2_i;
tilt_m = tilt_m_i;
z_s1 = z_s1_i;
z_s2 = z_s2_i;

n=1;

% Create random change matrix to module positions:
r_change = (r_var_max-r_var_min).*rand(runs,12) + r_var_min;
tilt_change = (tilt_var_max-tilt_var_min).*rand(runs,12) + tilt_var_min;
z_change = (z_var_max-z_var_min).*rand(runs,12) + z_var_min;
 
for g = 1:runs
%while n <= 12
%%%% Initialize loop parameters:
    
    % Loop:
    for k = 1:12
        %if k == 1
        %    r_s1(k+1) = r_s1_i(k+1) - r_tol(2,k);
        %    r_s2(k+1) = r_s2_i(k+1) - r_tol(2,k);
        %    z_s1(k+1) = z_s1_i(k+1) + z_tol(1,k);
        %    z_s2(k+1) = z_s2_i(k+1) + z_tol(1,k);
        %elseif k > 1
        %    r_s1(k+1) = r_s1_i(k+1) - r_tol(2,k);
        %    r_s2(k+1) = r_s2_i(k+1) - r_tol(2,k);
        %    z_s1(k+1) = z_s1_i(k+1) + z_tol(1,k);
        %    z_s2(k+1) = z_s2_i(k+1) + z_tol(1,k);

        %    r_s1(k) = r_s1_i(k) + r_tol(1,k-1);
        %    r_s2(k) = r_s2_i(k) + r_tol(1,k-1);
        %    z_s1(k) = z_s1_i(k) - z_tol(2,k-1);
        %    z_s2(k) = z_s2_i(k) - z_tol(2,k-1);
        %end
        
        % Calculate z & y of points 1, 2, 3 & 4

        % y-coordinates:
        p1_y = r_s1(k+1) + ( (h_m / 2) * sin(tilt_m(k+1)) );
        p2_y = r_s2(k+1) + ( (h_m / 2) * sin(tilt_m(k+1)) );
        p3_y = r_s1(k) - ( (h_m / 2) * sin(tilt_m(k)) );
        p4_y = r_s2(k) - ( (h_m / 2) * sin(tilt_m(k)) );
        % z-coordinates:
        p1_z = z_s1(k+1) - ( (h_m / 2) * cos(tilt_m(k+1)) );
        p2_z = z_s2(k+1) - ( (h_m / 2) * cos(tilt_m(k+1)) );
        p3_z = z_s1(k) + ( (h_m / 2) * cos(tilt_m(k)) );
        p4_z = z_s2(k) + ( (h_m / 2) * cos(tilt_m(k)) );

        % Fit polynomials of the first order to point 1-3, 1-4, 2-3 and 2-4:
        P_13 = polyfit([p1_z, p3_z], [p1_y, p3_y], 1); % 1-3
        P_14 = polyfit([p1_z, p4_z], [p1_y, p4_y], 1); % 1-4
        P_23 = polyfit([p2_z, p3_z], [p2_y, p3_y], 1); % 2-3
        P_24 = polyfit([p2_z, p4_z], [p2_y, p4_y], 1); % 2-4

        % Calculate z & y of points 5, 6, 7 & 8

        % y-coordinates:
        p5_y = r_s1(k+1) - ( (h_m / 2) * sin(tilt_m(k+1)) );
        p6_y = r_s2(k+1) - ( (h_m / 2) * sin(tilt_m(k+1)) );
        p7_y = r_s1(k) + ( (h_m / 2) * sin(tilt_m(k)) );
        p8_y = r_s2(k) + ( (h_m / 2) * sin(tilt_m(k)) );
        % z-coordinates:
        p5_z = z_s1(k+1) + ( (h_m / 2) * cos(tilt_m(k+1)) );
        p6_z = z_s2(k+1) + ( (h_m / 2) * cos(tilt_m(k+1)) );
        p7_z = z_s1(k) - ( (h_m / 2) * cos(tilt_m(k)) );
        p8_z = z_s2(k) - ( (h_m / 2) * cos(tilt_m(k)) );

        % Fit polynomials of the first order to point 1-5, 2-6, 3-7 and 4-8:
        P_15 = polyfit([p1_z, p5_z], [p1_y, p5_y], 1); % 1-5
        P_26 = polyfit([p2_z, p6_z], [p2_y, p6_y], 1); % 2-6
        P_37 = polyfit([p3_z, p7_z], [p3_y, p7_y], 1); % 3-7
        P_48 = polyfit([p4_z, p8_z], [p4_y, p8_y], 1); % 4-8

        % Calculate intersections of track polynomials and sensor polynomials:
        coe_resultants_P13_1 = P_13 - P_26; % sensor 2 in module i+1
        coe_resultants_P13_2 = P_13 - P_48; % sensor 1 in module i
        coe_resultants_P14_1 = P_14 - P_26; % sensor 2 in module i+1
        coe_resultants_P14_2 = P_14 - P_37; % sensor 1 in module i
        coe_resultants_P23_1 = P_23 - P_15; % sensor 1 in module i+1
        coe_resultants_P23_2 = P_23 - P_48; % sensor 2 in module i
        coe_resultants_P24_1 = P_24 - P_15; % sensor 1 in module i+1
        coe_resultants_P24_2 = P_24 - P_37; % sensor 1 in module i

        root_P13_1 = roots(coe_resultants_P13_1);
        root_P13_2 = roots(coe_resultants_P13_2);
        root_P14_1 = roots(coe_resultants_P14_1);
        root_P14_2 = roots(coe_resultants_P14_2);
        root_P23_1 = roots(coe_resultants_P23_1);
        root_P23_2 = roots(coe_resultants_P23_2);
        root_P24_1 = roots(coe_resultants_P24_1);
        root_P24_2 = roots(coe_resultants_P24_2);
        fprintf('z-coverage %d: \n', k)
        
        % Z-overlap:
        % Fit polynomials from 0 to points 3 & 4:
        P_03 = polyfit([0, p3_z], [0, p3_y], 1); % 0-3
        P_04 = polyfit([0, p4_z], [0, p4_y], 1); % 0-4
        % Check where they intersect sensor 1 of k+1 module:
        overlap_P03_z = (P_15(2) - P_03(2)) / (P_03(1) - P_15(1)) 
        overlap_P03_y = polyval(P_03, overlap_P03_z);
        overlap_P04_z = (P_15(2) - P_04(2)) / (P_04(1) - P_15(1))
        overlap_P04_y = polyval(P_04, overlap_P04_z);
        overlap_P03 = sqrt(((overlap_P03_z-p1_z)^2) + (p1_y-overlap_P03_y)^2)
        overlap_P04 = sqrt(((overlap_P04_z-p1_z)^2) + (p1_y-overlap_P04_y)^2)
        
        % Check which polynomial intersects all sensors:
        if (root_P13_1(1) >= p2_z) && (root_P13_2(1) <= p4_z)
            % Mark z-coverage as the polynomials z-value at y=0:
            z_cov(k) = (-P_13(2)) / P_13(1);
            z_overlap(k) = overlap_P03;
            %fprintf('Polynomial used: P13 \n\n')
        end

        if (root_P14_1(1) >= p2_z) && (root_P14_2(1) <= p3_z)
            z_cov(k) = (-P_14(2)) / P_14(1);
            z_overlap(k) = overlap_P04;
            %fprintf('Polynomial used: P14 \n\n')
        end

        if (root_P23_1(1) >= p1_z) && (root_P23_2(1) <= p4_z)
            z_cov(k) = (-P_23(2)) / P_23(1);
            z_overlap(k) = overlap_P03;
            %fprintf('Polynomial used: P23 \n\n')
        end

        if (root_P24_1(1) >= p1_z) && (root_P24_2(1) <= p3_z)
            z_cov(k) = (-P_24(2)) / P_24(1);
            z_overlap(k) = overlap_P04;
            %fprintf('Polynomial used: P24 \n\n')
        end
        fprintf('z-coverage with P13: %d\n', (-P_13(2)) / P_13(1))
        fprintf('z-coverage with P14: %d\n', (-P_14(2)) / P_14(1))
        fprintf('z-coverage with P23: %d\n', (-P_23(2)) / P_23(1))
        fprintf('z-coverage with P24: %d\n', (-P_24(2)) / P_24(1))
    end
    
    % Save results:
    %z_cov_results(g,1:12) = z_cov;
    % Displace modules incrementally:
end