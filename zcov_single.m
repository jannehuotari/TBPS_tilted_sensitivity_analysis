function [z_cov_single_result] = zcov_single(layer, modules, index, z_change, r_change, tilt_change, runs_n)
%% Geometry initialization
% Values from tkLayout with current geometry

h_m = 46.944; % Height of module

%%%%%% Layer 1 %%%%%%

% Layer 1 module sensor spacing:
sensor_spacing_layer1 = [2.6, 2.6, 2.6, 2.6, 4, 4, 4, 4, 4, 4, 4, 4, 4];
% Tilt angle in radians:
tilt_layer1 = [0, 47, 47, 47, 60, 60, 60, 60, 74, 74, 74, 74, 74] .* (pi/180);

%%% Layer 1 inner modules:

% Inner module radius:
r_m_layer1_inner = [217.537, 252, 252, 252, 249, 249, 249, 249, 248, 248, 248, 248, 248];
% Inner module position in z:
z_m_layer1_inner = [128.763, 172.095, 217.118, 267.987, 315.184, 374.203, 443.696, 526.772, 611.087, 720.246, 850.617, 1004.205, 1182.332];

% Inner module sensor positions in z:
z_s1_layer1_inner = zeros(1,13); 
z_s2_layer1_inner = zeros(1,13);
for n = 1:13
    z_s1_layer1_inner(n) = z_m_layer1_inner(n) + ( (sensor_spacing_layer1(n) * sin(tilt_layer1(n))) / 2 );
    z_s2_layer1_inner(n) = z_m_layer1_inner(n) - ( (sensor_spacing_layer1(n) * sin(tilt_layer1(n))) / 2 );
end

% Inner module sensor positions in r:
r_s1_layer1_inner = zeros(1,13);
r_s2_layer1_inner = zeros(1,13);
for m = 1:13
    r_s1_layer1_inner(m) = r_m_layer1_inner(m) + ( (sensor_spacing_layer1(m) * cos(tilt_layer1(m))) / 2 );
    r_s2_layer1_inner(m) = r_m_layer1_inner(m) - ( (sensor_spacing_layer1(m) * cos(tilt_layer1(m))) / 2 );
end

%%% Layer 1 outer modules:

% Outer module radius:
r_m_layer1_outer = [241.338 265, 265, 265, 259, 259, 259, 259, 254, 254, 254, 254, 254];
% Outer module position in z:
z_m_layer1_outer = [128.763, 181.891, 226.914, 277.784, 327.101, 386.121, 455.614, 538.690, 623.389, 732.548, 862.919, 1016.507, 1194.634];

% Outer module sensor positions in z:
z_s1_layer1_outer = zeros(1,13); 
z_s2_layer1_outer = zeros(1,13);
for n = 1:13
    z_s1_layer1_outer(n) = z_m_layer1_outer(n) + ( (sensor_spacing_layer1(n) * sin(tilt_layer1(n))) / 2 );
    z_s2_layer1_outer(n) = z_m_layer1_outer(n) - ( (sensor_spacing_layer1(n) * sin(tilt_layer1(n))) / 2 );
end

% Outer module sensor positions in r:
r_s1_layer1_outer = zeros(1,13);
r_s2_layer1_outer = zeros(1,13);
for m = 1:13
    r_s1_layer1_outer(m) = r_m_layer1_outer(m) + ( (sensor_spacing_layer1(m) * cos(tilt_layer1(m))) / 2 );
    r_s2_layer1_outer(m) = r_m_layer1_outer(m) - ( (sensor_spacing_layer1(m) * cos(tilt_layer1(m))) / 2 );
end

%%%%%% Layer 2 %%%%%%

% Layer 2 module sensor spacing:
sensor_spacing_layer2 = [1.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 4, 4, 4, 4, 4];
% Tilt angle in radians:
tilt_layer2 = [0, 40, 40, 40, 55, 55, 55, 55, 68, 68, 68, 68, 68] .* (pi/180);

%%% Layer 2 inner modules:

% Inner module radius:
r_m_layer2_inner = [346.330, 374, 374, 374, 372, 372, 372, 372, 371, 371, 371, 371, 371];
% Inner module position in z:
z_m_layer2_inner = [222.566, 269.888, 321.982, 378.406, 431.314, 495.745, 567.262, 646.108, 725.236, 820.642, 927.596, 1048.134, 1181.542];

% Inner module sensor positions in z:
z_s1_layer2_inner = zeros(1,13); 
z_s2_layer2_inner = zeros(1,13);
for n = 1:13
    z_s1_layer2_inner(n) = z_m_layer2_inner(n) + ( (sensor_spacing_layer2(n) * sin(tilt_layer2(n))) / 2 );
    z_s2_layer2_inner(n) = z_m_layer2_inner(n) - ( (sensor_spacing_layer2(n) * sin(tilt_layer2(n))) / 2 );
end

% Inner module sensor positions in r:
r_s1_layer2_inner = zeros(1,13);
r_s2_layer2_inner = zeros(1,13);
for m = 1:13
    r_s1_layer2_inner(m) = r_m_layer2_inner(m) + ( (sensor_spacing_layer2(m) * cos(tilt_layer2(m))) / 2 );
    r_s2_layer2_inner(m) = r_m_layer2_inner(m) - ( (sensor_spacing_layer2(m) * cos(tilt_layer2(m))) / 2 );
end

%%% Layer 2 outer modules:

% Outer module radius:
r_m_layer2_outer = [370.130, 387.000, 387.000, 387.000, 381.000, 381.000, 381.000, 381.000, 378.000, 378.000, 378.000, 378.000, 378.000];
% Outer module position in z:
z_m_layer2_outer = [222.566, 277.394, 329.488, 385.911, 440.314, 504.745, 576.262, 655.108, 736.439, 831.844, 938.798, 1059.336, 1192.744];

% Outer module sensor positions in z:
z_s1_layer2_outer = zeros(1,13); 
z_s2_layer2_outer = zeros(1,13);
for n = 1:13
    z_s1_layer2_outer(n) = z_m_layer2_outer(n) + ( (sensor_spacing_layer2(n) * sin(tilt_layer2(n))) / 2 );
    z_s2_layer2_outer(n) = z_m_layer2_outer(n) - ( (sensor_spacing_layer2(n) * sin(tilt_layer2(n))) / 2 );
end

% Outer module sensor positions in r:
r_s1_layer2_outer = zeros(1,13);
r_s2_layer2_outer = zeros(1,13);
for m = 1:13
    r_s1_layer2_outer(m) = r_m_layer2_outer(m) + ( (sensor_spacing_layer2(m) * cos(tilt_layer2(m))) / 2 );
    r_s2_layer2_outer(m) = r_m_layer2_outer(m) - ( (sensor_spacing_layer2(m) * cos(tilt_layer2(m))) / 2 );
end

%%%%%% Layer 3 %%%%%%

% Layer 3 module sensor spacing:
sensor_spacing_layer3 = [1.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6];
% Tilt angle in radians:
tilt_layer3 = [0, 44, 44, 44, 44, 44, 44, 60, 60, 60, 60, 60, 60] .* (pi/180);

%%% Layer 3 inner modules:

% Inner module radius:
r_m_layer3_inner = [498.630, 524.500, 524.500, 524.500, 524.500, 524.500, 524.500, 522.500, 522.500, 522.500, 522.500, 522.500, 522.500];
% Inner module position in z:
z_m_layer3_inner = [315.734, 363.047, 416.152, 472.684, 533.009, 597.236, 665.109, 727.825, 804.833, 888.179, 978.285, 1075.169, 1179.356];

% Inner module sensor positions in z:
z_s1_layer3_inner = zeros(1,13); 
z_s2_layer3_inner = zeros(1,13);
for n = 1:13
    z_s1_layer3_inner(n) = z_m_layer3_inner(n) + ( (sensor_spacing_layer3(n) * sin(tilt_layer3(n))) / 2 );
    z_s2_layer3_inner(n) = z_m_layer3_inner(n) - ( (sensor_spacing_layer3(n) * sin(tilt_layer3(n))) / 2 );
end

% Inner module sensor positions in r:
r_s1_layer3_inner = zeros(1,13);
r_s2_layer3_inner = zeros(1,13);
for m = 1:13
    r_s1_layer3_inner(m) = r_m_layer3_inner(m) + ( (sensor_spacing_layer3(m) * cos(tilt_layer3(m))) / 2 );
    r_s2_layer3_inner(m) = r_m_layer3_inner(m) - ( (sensor_spacing_layer3(m) * cos(tilt_layer3(m))) / 2 );
end

%%% Layer 3 outer modules:

% Outer module radius:
r_m_layer3_outer = [522.430, 536.500, 536.500, 536.500, 536.500, 536.500, 536.500, 530.500, 530.500, 530.500, 530.500, 530.500, 530.500];
% Outer module position in z:
z_m_layer3_outer = [315.734, 371.141, 424.246, 480.778, 541.103, 605.330, 673.203, 737.359, 814.367, 897.713, 987.819, 1084.703, 1188.890];

% Outer module sensor positions in z:
z_s1_layer3_outer = zeros(1,13); 
z_s2_layer3_outer = zeros(1,13);
for n = 1:13
    z_s1_layer3_outer(n) = z_m_layer3_outer(n) + ( (sensor_spacing_layer3(n) * sin(tilt_layer3(n))) / 2 );
    z_s2_layer3_outer(n) = z_m_layer3_outer(n) - ( (sensor_spacing_layer3(n) * sin(tilt_layer3(n))) / 2 );
end

% Outer module sensor positions in r:
r_s1_layer3_outer = zeros(1,13);
r_s2_layer3_outer = zeros(1,13);
for m = 1:13
    r_s1_layer3_outer(m) = r_m_layer3_outer(m) + ( (sensor_spacing_layer3(m) * cos(tilt_layer3(m))) / 2 );
    r_s2_layer3_outer(m) = r_m_layer3_outer(m) - ( (sensor_spacing_layer3(m) * cos(tilt_layer3(m))) / 2 );
end
%% Actual computation
runs = runs_n; % Simulation runs

% Initialize some arrays
z_cov_single_result = zeros(runs,1);

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

for g = 1:runs
    for t = 1:13
        if t == 1
            r_s1(t) = r_s1_i(t);
            r_s2(t) = r_s2_i(t);
            tilt_m(t) = tilt_m_i(t);
            z_s1(t) = z_s1_i(t);
            z_s2(t) = z_s2_i(t);
        elseif t > 1
            r_s1(t) = r_s1_i(t) + r_change(g,t-1);
            r_s2(t) = r_s2_i(t) + r_change(g,t-1);
            tilt_m(t) = tilt_m_i(t) + tilt_change(g,t-1);
            z_s1(t) = z_s1_i(t) + z_change(g,t-1);
            z_s2(t) = z_s2_i(t) + z_change(g,t-1);
        end
    end
        
        % Calculate z & y of points 1, 2, 3 & 4

        % y-coordinates:
        p1_y = r_s1(index+1) + ( (h_m / 2) * sin(tilt_m(index+1)) );
        p2_y = r_s2(index+1) + ( (h_m / 2) * sin(tilt_m(index+1)) );
        p3_y = r_s1(index) - ( (h_m / 2) * sin(tilt_m(index)) );
        p4_y = r_s2(index) - ( (h_m / 2) * sin(tilt_m(index)) );
        % z-coordinates:
        p1_z = z_s1(index+1) - ( (h_m / 2) * cos(tilt_m(index+1)) );
        p2_z = z_s2(index+1) - ( (h_m / 2) * cos(tilt_m(index+1)) );
        p3_z = z_s1(index) + ( (h_m / 2) * cos(tilt_m(index)) );
        p4_z = z_s2(index) + ( (h_m / 2) * cos(tilt_m(index)) );

        % Fit polynomials of the first order to point 1-3, 1-4, 2-3 and 2-4:
        P_13 = polyfit([p1_z, p3_z], [p1_y, p3_y], 1); % 1-3
        P_14 = polyfit([p1_z, p4_z], [p1_y, p4_y], 1); % 1-4
        P_23 = polyfit([p2_z, p3_z], [p2_y, p3_y], 1); % 2-3
        P_24 = polyfit([p2_z, p4_z], [p2_y, p4_y], 1); % 2-4

        % Calculate z & y of points 5, 6, 7 & 8

        % y-coordinates:
        p5_y = r_s1(index+1) - ( (h_m / 2) * sin(tilt_m(index+1)) );
        p6_y = r_s2(index+1) - ( (h_m / 2) * sin(tilt_m(index+1)) );
        p7_y = r_s1(index) + ( (h_m / 2) * sin(tilt_m(index)) );
        p8_y = r_s2(index) + ( (h_m / 2) * sin(tilt_m(index)) );
        % z-coordinates:
        p5_z = z_s1(index+1) + ( (h_m / 2) * cos(tilt_m(index+1)) );
        p6_z = z_s2(index+1) + ( (h_m / 2) * cos(tilt_m(index+1)) );
        p7_z = z_s1(index) - ( (h_m / 2) * cos(tilt_m(index)) );
        p8_z = z_s2(index) - ( (h_m / 2) * cos(tilt_m(index)) );

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
        %fprintf('z-coverage %d: \n', k)
        % Check which polynomial intersects all sensors:
        if (root_P13_1(1) >= p2_z) && (root_P13_2(1) <= p4_z)
            % Mark z-coverage as the polynomials z-value at y=0:
            z_cov = (-P_13(2)) / P_13(1);
            %fprintf('Polynomial used: P13 \n\n')
        end

        if (root_P14_1(1) >= p2_z) && (root_P14_2(1) <= p3_z)
            z_cov = (-P_14(2)) / P_14(1);
            %fprintf('Polynomial used: P14 \n\n')
        end

        if (root_P23_1(1) >= p1_z) && (root_P23_2(1) <= p4_z)
            z_cov = (-P_23(2)) / P_23(1);
            %fprintf('Polynomial used: P23 \n\n')
        end

        if (root_P24_1(1) >= p1_z) && (root_P24_2(1) <= p3_z)
            z_cov = (-P_24(2)) / P_24(1);
            %fprintf('Polynomial used: P24 \n\n')
        end
        %fprintf('z-coverage with P13: %d\n', (-P_13(2)) / P_13(1))
        %fprintf('z-coverage with P14: %d\n', (-P_14(2)) / P_14(1))
        %fprintf('z-coverage with P23: %d\n', (-P_23(2)) / P_23(1))
        %fprintf('z-coverage with P24: %d\n', (-P_24(2)) / P_24(1))
    
    % Save results:
    z_cov_single_result(g,1) = z_cov;
end

end