% function returns phi-overlap value for the specified layer and ring.
% layer = layer to analyze (integer value).

% ring = ring to analyze (integer value).

% z_change, r_change etc... = runs x 12 array of changes applied. Values on
% certain array columns apply to that specific ring. For example, value in
% column 5 of a change array will be used to modify ring 5 in specified
% layer.

% z_overlap = runs x 12 array of z_overlap values associated with the
% modules. Notice that this is the actual z-overlap value, rather than a
% change from nominal z-overlap value.

% runs_n = integet of how many times results are ran. Same size as the
% amount of rows in change arrays.

function [phioverlap] = phi_overlap(layer, ring, z_change, r_change, x_change, x_tilt_change, y_tilt_change, z_tilt_change, z_overlap, runs_n)

%% Testing:
%clear;
%layer = 1;
%ring = 1;
%z_change = 0;
%r_change = 0;
%x_change = 0;
%x_tilt_change = 0;
%y_tilt_change = 0;
%z_tilt_change = 0;
%z_overlap = 10;
%runs_n = 1;
%%
% Phi-overlap calculation

%% Geometry initialization
% Values from tkLayout with current geometry

h_m = 46.944; % Height of module
w_m = 96.000; % Width of module
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

% Unique code starts here

phioverlap = zeros(runs_n,1);

for k = 1:runs_n
    % Initialize geometrical values and apply changes to outer module (errors)
    if layer == 1
        r_inner = r_s1_layer1_inner(ring+1);
        z_inner = z_s1_layer1_inner(ring+1);
        r_outer = r_s1_layer1_outer(ring+1)+r_change(k, ring);
        z_outer = z_s1_layer1_outer(ring+1)+z_change(k, ring);
        spacing = sensor_spacing_layer1(ring+1);
        tilt_x = tilt_layer1(ring+1);
        tilt_y = 0;
        phi = 20*(pi/180);
    elseif layer == 2
        r_inner = r_s1_layer2_inner(ring+1);
        z_inner = z_s1_layer2_inner(ring+1);
        r_outer = r_s1_layer2_outer(ring+1)+r_change(k, ring);
        z_outer = z_s1_layer2_outer(ring+1)+z_change(k, ring);
        spacing = sensor_spacing_layer2(ring+1);
        tilt_x = tilt_layer2(ring+1);
        tilt_y = 0;
        phi = 13.846154*(pi/180);
    elseif layer == 3
        r_inner = r_s1_layer3_inner(ring+1);
        z_inner = z_s1_layer3_inner(ring+1);
        r_outer = r_s1_layer3_outer(ring+1)+r_change(k, ring);
        z_outer = z_s1_layer3_outer(ring+1)+z_change(k, ring);
        spacing = sensor_spacing_layer3(ring+1);
        tilt_x = tilt_layer3(ring+1);
        tilt_y = 0;
        phi = 10*(pi/180);
    end
    
    % Calculate phi-overlap:
    
    % Refer to module sensitivity analysis report in sharepoint for further
    % clarification (phi-overlap section and appendix B).
    
    % Figure out 2d coordinates of point 6 in module local coordinate
    % system:
    % Position in local module coordinate system with phi = 0:
    
    % Info on parameter naming:
    % l = module local coordinate system
    % l_t = module local coordinate system with module rotated to right
    % position
    % real = z-overlap dimension not taken into account
    
    p6_x_l = -(w_m/2) * cos(tilt_y+y_tilt_change(k,ring));
    p6_y_l = ((h_m/2)-z_overlap(k, ring)) * sin(tilt_x+x_tilt_change(k, ring));
    
    p6_real_x_l = -(w_m/2) * cos(tilt_y+y_tilt_change(k,ring));
    p6_real_y_l = h_m/2 * sin(tilt_x+x_tilt_change(k, ring));
    
    p5_real_x_l = (w_m/2) * cos(tilt_y+y_tilt_change(k,ring));
    p5_real_y_l = h_m/2 * sin(tilt_x+x_tilt_change(k, ring));
    
    p7_real_x_l = -(w_m/2) * cos(tilt_y+y_tilt_change(k,ring));
    p7_real_y_l = -h_m/2 * sin(tilt_x+x_tilt_change(k, ring));
    
    % Position in local module coordinate system with module rotated according to phi:
    
    p6_x_l_t = p6_x_l*cos(phi+z_tilt_change(k, ring)) + p6_y_l*sin(phi+z_tilt_change(k, ring));
    p6_y_l_t = -p6_x_l*sin(phi+z_tilt_change(k, ring)) + p6_y_l*cos(phi+z_tilt_change(k, ring));
    
    p6_real_x_l_t = p6_real_x_l*cos(phi+z_tilt_change(k, ring)) + p6_real_y_l*sin(phi+z_tilt_change(k, ring));
    p6_real_y_l_t = -p6_real_x_l*sin(phi+z_tilt_change(k, ring)) + p6_real_y_l*cos(phi+z_tilt_change(k, ring));
    
    p5_real_x_l_t = p5_real_x_l*cos(phi+z_tilt_change(k, ring)) + p5_real_y_l*sin(phi+z_tilt_change(k, ring));
    p5_real_y_l_t = -p5_real_x_l*sin(phi+z_tilt_change(k, ring)) + p5_real_y_l*cos(phi+z_tilt_change(k, ring));
    
    p7_real_x_l_t = p7_real_x_l*cos(phi+z_tilt_change(k, ring)) + p7_real_y_l*sin(phi+z_tilt_change(k, ring));
    p7_real_y_l_t = -p7_real_x_l*sin(phi+z_tilt_change(k, ring)) + p7_real_y_l*cos(phi+z_tilt_change(k, ring));
    
    % Figure out point coordinates in ring coordinate system:
    
    % Coordinates of points in ring coordinate system:
    p1_x_s1 = (w_m/2) * cos(tilt_y);
    p1_y_s1 = r_inner + ((h_m/2)-z_overlap(k, ring)) * sin(tilt_x);
    p1_z_s1 = z_inner - ((h_m/2)-z_overlap(k, ring)) * cos(tilt_x);
   
    p2_x_s1 = -(w_m/2) * cos(tilt_y);
    p2_y_s1 = r_inner + ((h_m/2)-z_overlap(k, ring)) * sin(tilt_x);
    p2_z_s1 = z_inner - ((h_m/2)-z_overlap(k, ring)) * cos(tilt_x);
    
    p4_x_s1 = (w_m/2) * cos(tilt_y);
    p4_y_s1 = r_inner - (h_m/2) * sin(tilt_x);
    p4_z_s1 = z_inner + (h_m/2) * cos(tilt_x);
    
    p6_x_s1 = r_outer * cos((pi/2) - phi) + p6_x_l_t + x_change(k,ring)*sin(phi+z_tilt_change(k, ring));
    p6_y_s1 = r_outer * sin((pi/2) - phi) + p6_y_l_t - x_change(k,ring)*cos(phi+z_tilt_change(k, ring));
    p6_z_s1 = z_outer - ((h_m/2)-z_overlap(k, ring)) * cos(tilt_x + x_tilt_change(k, ring));
    
    p6_real_x_s1 = r_outer * cos((pi/2) - phi) + p6_real_x_l_t + x_change(k,ring)*sin(phi+z_tilt_change(k, ring));
    p6_real_y_s1 = r_outer * sin((pi/2) - phi) + p6_real_y_l_t - x_change(k,ring)*cos(phi+z_tilt_change(k, ring));
    p6_real_z_s1 = z_outer - (h_m/2) * cos(tilt_x + x_tilt_change(k, ring));
    
    p5_real_x_s1 = r_outer * cos((pi/2) - phi) + p5_real_x_l_t + x_change(k,ring)*sin(phi+z_tilt_change(k, ring));
    p5_real_y_s1 = r_outer * sin((pi/2) - phi) + p5_real_y_l_t - x_change(k,ring)*cos(phi+z_tilt_change(k, ring));
    p5_real_z_s1 = z_outer - (h_m/2) * cos(tilt_x + x_tilt_change(k, ring));
    
    p7_real_x_s1 = r_outer * cos((pi/2) - phi) + p7_real_x_l_t + x_change(k,ring)*sin(phi+z_tilt_change(k, ring));
    p7_real_y_s1 = r_outer * sin((pi/2) - phi) + p7_real_y_l_t - x_change(k,ring)*cos(phi+z_tilt_change(k, ring));
    p7_real_z_s1 = z_outer + (h_m/2) * cos(tilt_x + x_tilt_change(k, ring));
    
    % Calculation for phi-overlap that is determined from a particle that
    % traverses from origin to point 6:
    % Vector from origo to point 6:
    op6 = [p6_x_s1 p6_y_s1 p6_z_s1];
    % Vector from p1 to p2
    p1p2 = [p2_x_s1-p1_x_s1 p2_y_s1-p1_y_s1 p2_z_s1-p1_z_s1];
    % Vector prom p1 to p4
    p1p4 = [p4_x_s1-p1_x_s1 p4_y_s1-p1_y_s1 p4_z_s1-p1_z_s1];
    % Vector normal to p1p2 and p1p4 (normal of plane):
    n_plane_left_s = cross(p1p2,p1p4);
    % Figure out at which t, op6 intersects plane of left sensor:
    t_op6 = (n_plane_left_s(1)*p1_x_s1 + n_plane_left_s(2)*p1_y_s1 + n_plane_left_s(3)*p1_z_s1) / (n_plane_left_s(1)*op6(1) + n_plane_left_s(2)*op6(2) + n_plane_left_s(3)*op6(3));
    
    intersection_op6 = op6.*t_op6;
    
    % Calculation for phi-overlap that is determined from a particle that
    % traverses from origin to point 1:
    op1 = [p1_x_s1 p1_y_s1 p1_z_s1];
    % Vector from p6 to p5
    p6p5 = [p5_real_x_s1-p6_real_x_s1 p5_real_y_s1-p6_real_y_s1 p5_real_z_s1-p6_real_z_s1];
    % Vector prom p6 to p7
    p6p7 = [p7_real_x_s1-p6_real_x_s1 p7_real_y_s1-p6_real_y_s1 p7_real_z_s1-p6_real_z_s1];
    % Vector normal to p6p5 and p6p7 (normal of plane):
    n_plane_right_s = cross(p6p5,p6p7);
    % Figure out at which t, op1 intersects plane of right sensor:
    t_op1 = (n_plane_right_s(1)*p6_real_x_s1 + n_plane_right_s(2)*p6_real_y_s1 + n_plane_right_s(3)*p6_real_z_s1) / (n_plane_right_s(1)*op1(1) + n_plane_right_s(2)*op1(2) + n_plane_right_s(3)*op1(3));
    
    intersection_op1 = op1.*t_op1;
    
    % Calculate phi-overlaps from intersections; calculated as intersection
    % points distance from sensors edge in x for right sensor, for left sensor, just the distance from x coordinate of point 1:
    
    % Pay no attention to this:
    % Formulate left sensors edge as a line from point 1 to point 4:
    %left_sensor_edge = polyfit([p1_x_s1, p4_x_s1], [p1_y_s1, p4_y_s1], 1);
    
    % Formulate right sensors edge as a line from point 6 to point 7:
    right_sensor_edge = polyfit([p6_real_x_s1, p7_real_x_s1], [p6_real_y_s1, p7_real_y_s1], 1);
    
    % Left sensor edges x coordinate at op6 intersection point y:
    %left_sensor_edge_x = (intersection_op6(2) - left_sensor_edge(2)) / left_sensor_edge(1);
    % Right sensor edges x coordinate at op1 intersection point y:
    right_sensor_edge_x = (intersection_op1(2) - right_sensor_edge(2)) / right_sensor_edge(1);
    
    %phi-overlaps:
    phi_overlap_op6 = p1_x_s1 - intersection_op6(1);
    phi_overlap_op1 = intersection_op1(1) - right_sensor_edge_x;
    % Check with which trajectory (op1 or op6) phi-overlap is smaller (actual phi-overlap will be the one smaller):
    if runs_n > 1
        if phi_overlap_op6 < phi_overlap_op1
            phioverlap(k, 1) = phi_overlap_op6;
        elseif phi_overlap_op6 >= phi_overlap_op1
            phioverlap(k, 1) = phi_overlap_op1;
        end
    elseif runs_n == 1
        if phi_overlap_op6 < phi_overlap_op1
            phioverlap = phi_overlap_op6;
        elseif phi_overlap_op6 >= phi_overlap_op1
            phioverlap = phi_overlap_op1;
        end
    end
    
end
end