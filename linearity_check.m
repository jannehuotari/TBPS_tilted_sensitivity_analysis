% The purpose of this .m file is to illustrate the relationship between the
% dimensional parameters and z-coverage / phi-overlap. Most of the
% relationships are linear, except for alpha variation.

%% Linearity check
% Initialization
clear;

z_change_z = repmat(linspace(0,5,100)', 1, 12);
r_change_z = repmat(linspace(0,5,100)', 1, 12);
tilt_change_z = repmat(linspace(0,deg2rad(1),100)', 1, 12);
null_z = zeros(1, 100);
z_change_p = repmat(linspace(-0.4,0.4,100)', 1, 12);
r_change_p = repmat(linspace(0,0.8,100)', 1, 12);
x_change_p = repmat(linspace(0,0.8,100)', 1, 12);
x_tilt_change_p = repmat(linspace(0,deg2rad(1),100)', 1, 12);
y_tilt_change_p = repmat(linspace(0,deg2rad(1),100)', 1, 12);
z_tilt_change_p = repmat(linspace(0,deg2rad(1),100)', 1, 12);
z_overlap_p = repmat(linspace(0,5,100)', 1, 12);
null_p = zeros(100, 12);

z_result = zeros(1, 100);
r_z = zeros(1, 100);
tilt_z = zeros(1, 100);

z_p = zeros(1, 100);
r_p = zeros(1, 100);
x_p = zeros(1, 100);
x_tilt_p = zeros(1, 100);
y_tilt_p = zeros(1, 100);
z_tilt_p = zeros(1, 100);
z_ove_p = zeros(1, 100);

%zcov(layer, modules, z_change, r_change, tilt_change, runs_n)

% z_change
for i = 1:100
    temp = zcov(1, 'inner', z_change_z(i, 1:12), null_z, null_z, 1);
    z_result(i) = temp(1);
end
figure;
grid on;
plot(z_change_z, z_result);
xlabel('z-change');
ylabel('z-coverage');
title('z-change-z-coverage relationship');
% r_change
for i = 1:100
    temp = zcov(1, 'inner', null_z, r_change_z(i,1:12), null_z, 1);
    z_result(i) = temp(1);
end
figure;
grid on;
plot(r_change_z, z_result);
xlabel('r-change');
ylabel('z-coverage');
title('r-change-z-coverage relationship');
%tilt_change
for i = 1:100
    temp = zcov(1, 'inner', null_z, null_z, tilt_change_z(i,1:12), 1);
    z_result(i) = temp(1);
end
figure;
grid on;
plot(tilt_change_z, z_result);
xlabel('tilt-change');
ylabel('z-coverage');
title('tilt-change-z-coverage relationship');

%phi_overlap(layer, ring, z_change, r_change, x_change, x_tilt_change, y_tilt_change, z_tilt_change, z_overlap, runs_n)

% z_change
for i = 1:100
    temp = phi_overlap(1, 1, z_change_p(i,1:12), null_p, null_p, null_p, null_p, null_p, null_p, 1);
    p_result(i) = temp;
end
figure;
grid on;
plot(z_change_p, p_result);
xlabel('z-change');
ylabel('phi-overlap');
title('z-change-phi-overlap relationship');
% r_change
for i = 1:100
    temp = phi_overlap(1, 1, null_p, r_change_p(i,1:12), null_p, null_p, null_p, null_p, null_p, 1);
    p_result(i) = temp(1);
end
figure;
grid on;
plot(r_change_p, p_result);
xlabel('r-change');
ylabel('phi-overlap');
title('r-change-phi-overlap relationship');
%v_change
for i = 1:100
    temp = phi_overlap(1, 1, null_p, null_p, x_change_p(i,1:12), null_p, null_p, null_p, null_p, 1);
    p_result(i) = temp(1);
end
figure;
grid on;
plot(x_change_p, p_result);
xlabel('v-change');
ylabel('phi-overlap');
title('v-change-phi-overlap relationship');
%x_tilt_change
for i = 1:100
    temp = phi_overlap(1, 1, null_p, null_p, null_p, x_tilt_change_p(i,1:12), null_p, null_p, null_p, 1);
    p_result(i) = temp(1);
end
figure;
grid on;
plot(x_tilt_change_p, p_result);
xlabel('beta-change');
ylabel('phi-overlap');
title('beta-change-phi-overlap relationship');
%y_tilt_change
for i = 1:100
    temp = phi_overlap(1, 1, null_p, null_p, null_p, null_p, y_tilt_change_p(i,1:12), null_p, null_p, 1);
    p_result(i) = temp(1);
end
figure;
grid on;
plot(y_tilt_change_p, p_result);
xlabel('alpha-change');
ylabel('phi-overlap');
title('alpha-change-phi-overlap relationship');
%z_tilt_change
for i = 1:100
    temp = phi_overlap(1, 1, null_p, null_p, null_p, null_p, null_p, z_tilt_change_p(i,1:12), null_p, 1);
    p_result(i) = temp(1);
end
figure;
grid on;
plot(z_tilt_change_p, p_result);
xlabel('gamma-change');
ylabel('phi-overlap');
title('gamma-change-phi-overlap relationship');
%z_overlap_change
for i = 1:100
    temp = phi_overlap(1, 1, null_p, null_p, null_p, null_p, null_p, null_p, z_overlap_p(i,1:12), 1);
    p_result(i) = temp(1);
end
figure;
grid on;
plot(z_overlap_p, p_result);
xlabel('z-overlap-change');
ylabel('phi-overlap');
title('z-overlap-change-phi-overlap relationship');