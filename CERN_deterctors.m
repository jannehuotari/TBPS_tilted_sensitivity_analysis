% This is the main .m file from which the simulations can be ran from. To
% run simulations, run the initialization part here in the beginning at
% first. You can select values based on what you wish to simulate. Each
% type of simulation is divided into segments with descriptive titles.

%% Initialization
clear;
layer   = 3; % Layer to analyze (This applies to certain simulations where layers are analyzed individually)
modules = 'inner'; % "inner" or "outer" modules to analyze (This applies to certain simulations where inner / outer modules are analyzed individually)
allowed_cov = 70; % Allowed min z-coverance
runs    = 10000; % The amount of runs has to do with sensitivity analysis sections
tilt_tol = repmat(deg2rad(1),2,12); % Assumed module tile angle changes applied in regards to z-coverage
zfilename = 'z_coverage_tolerances'; % Name of the excel file to which code writes z-coverage tolerances
pfilename = 'phi_overlap_tolerances'; % Name of the excel file to which code writes phi-overlap tolerances
phi_results = 'worst'; % Use 'nominal' or 'worst' -case z-overlap when calculating phi-overlap related tolerances
%% z-coverage Tolerances
% Running this will create an xml file to the current folder specified in
% Matlab (<---- on the left bar)
% Tolerance estimation
z_change = linspace(0,50,10000);
r_change = linspace(0,50,10000);
ratio = 0.6; % Ratio of z-tolerance to r-tolerance

% Create xls file base for tolerance values:
zcovxls(zfilename);

% Fill with tolerance values:
for i = 1:6
    if i == 1
        layer = 1;
        modules = 'inner';
        z_tol = zcov_ztol(layer, modules, z_change, allowed_cov);
        r_tol = zcov_rtol(layer, modules, z_tol, r_change, allowed_cov, ratio);
        z_tol = z_tol.*ratio;
        xlswrite(zfilename, z_tol, 'Layer 1', 'A4:L5');
        xlswrite(zfilename, r_tol, 'Layer 1', 'A7:L8');
    elseif i == 2
        layer = 1;
        modules = 'outer';
        z_tol = zcov_ztol(layer, modules, z_change, allowed_cov);
        r_tol = zcov_rtol(layer, modules, z_tol, r_change, allowed_cov, ratio);
        z_tol = z_tol.*ratio;
        xlswrite(zfilename, z_tol, 'Layer 1', 'A11:L12');
        xlswrite(zfilename, r_tol, 'Layer 1', 'A14:L15');
    elseif i == 3
        layer = 2;
        modules = 'inner';
        z_tol = zcov_ztol(layer, modules, z_change, allowed_cov);
        r_tol = zcov_rtol(layer, modules, z_tol, r_change, allowed_cov, ratio);
        z_tol = z_tol.*ratio;
        xlswrite(zfilename, z_tol, 'Layer 2', 'A4:L5');
        xlswrite(zfilename, r_tol, 'Layer 2', 'A7:L8');
    elseif i == 4
        layer = 2;
        modules = 'outer';
        z_tol = zcov_ztol(layer, modules, z_change, allowed_cov);
        r_tol = zcov_rtol(layer, modules, z_tol, r_change, allowed_cov, ratio);
        z_tol = z_tol.*ratio;
        xlswrite(zfilename, z_tol, 'Layer 2', 'A11:L12');
        xlswrite(zfilename, r_tol, 'Layer 2', 'A14:L15');
    elseif i == 5
        layer = 3;
        modules = 'inner';
        z_tol = zcov_ztol(layer, modules, z_change, allowed_cov);
        r_tol = zcov_rtol(layer, modules, z_tol, r_change, allowed_cov, ratio);
        z_tol = z_tol.*ratio;
        xlswrite(zfilename, z_tol, 'Layer 3', 'A4:L5');
        xlswrite(zfilename, r_tol, 'Layer 3', 'A7:L8');
    elseif i == 6
        layer = 3;
        modules = 'outer';
        z_tol = zcov_ztol(layer, modules, z_change, allowed_cov);
        r_tol = zcov_rtol(layer, modules, z_tol, r_change, allowed_cov, ratio);
        z_tol = z_tol.*ratio;
        xlswrite(zfilename, z_tol, 'Layer 3', 'A11:L12');
        xlswrite(zfilename, r_tol, 'Layer 3', 'A14:L15');
    end
end
%% Phi-overlap sensitivity analysis
% Initialize wost-case z-overlap, nominal phi-overlap and worst case phi-overlap arrays:
phi = zeros(runs,12);

% Calculate z-overlaps with nominal module positions:
% Initialize applied dimensional changes as zero:
z_change = zeros(1,12);
r_change = zeros(1,12);
tilt_change = zeros(1,12);
z_overlap_varied = zeros(runs,12);

z_cov_result = zcov(layer, 'outer', z_change, r_change, tilt_change, 1);
z_overlap = z_cov_result(13:24);

% Implement small variance in dimensions:
z_change        = (0.8.*lhs(runs,12)-0.4);
r_change        = (0.8.*lhs(runs,12)-0.4);
x_change        = (0.8.*lhs(runs,12)-0.4);
x_tilt_change   = (2*lhs(runs,12)-1).* pi/180;
y_tilt_change   = (2.*lhs(runs,12)-1).* pi/180;
z_tilt_change   = (2.*lhs(runs,12)-1).* pi/180;
% z_overlap_change dependent on amount of z-coverage in layer. To simulate
% this, change higher in layer 1 and 2 than in layer 3 (layer 3 z-coverage very close to 70)
if layer == 1
    z_overlap_change = (6.*lhs(runs,12)-3);
elseif layer == 2
    z_overlap_change = (4.*lhs(runs,12)-2);
elseif layer == 3
    z_overlap_change = (1.*lhs(runs,12)-0.5);
end

% Apply z_overlap_change to nominal z_overlaps:
for i = 1:runs
    for k = 1:12
        z_overlap_varied(i, k) = z_overlap(k) + z_overlap_change(i, k);
    end
end

% Loop phi-overlap with randomized module positional:
for index=1:12
    phi = phi_overlap(layer, index,  z_change, r_change, x_change, x_tilt_change, y_tilt_change, z_tilt_change, z_overlap_varied, runs);
    Y = phi;
    % See how change vectors correlate with result vector:
    St = plcc([z_change(:,index), r_change(:,index), x_change(:,index), x_tilt_change(:,index), y_tilt_change(:,index), z_tilt_change(:,index), z_overlap_change(:,index)],Y);
    z_st(index) = St(1);
    r_st(index) = St(2);
    x_st(index) = St(3);
    x_t_st(index) = St(4);
    y_t_st(index) = St(5);
    z_t_st(index) = St(6);
    z_o_st(index) = St(7);
end
% Plot:
figure;
plot(1:12,z_st); grid on; hold on;
plot(1:12,r_st); plot(1:12,x_st); plot(1:12,x_t_st); plot(1:12,y_t_st); plot(1:12,z_t_st); plot(1:12,z_o_st);
ylabel('Correlation coefficient'); xlabel('Ring');
legend('z variation', 'r variation', 'v variation', 'beta variation', 'alpha variation', 'gamma variation', 'z-overlap variation');
title(['phi-overlap tolerance sensitivies with fixed tolerances - Layer ', num2str(layer)]);

%% Phi-overlap Tolerances
% Running this will create an xml file to the current folder specified in
% Matlab (<---- on the left bar)
% Tolerance estimation
z_change = linspace(0,250,100000);
r_change = linspace(0,250,100000);
x_change = linspace(0,250,100000);
x_tilt_change = linspace(0,deg2rad(360),100000);
z_tilt_change = linspace(0,deg2rad(360),100000);
z_ratio = 0.7;
r_ratio = 0.7;
x_ratio = 0.5;
x_tilt_ratio = 0.5;
z_tilt_ratio = 0.5;
limit = 0;
tilt_change = zeros(1,12);
z_change_zo = zeros(1,12);
r_change_zo = zeros(1,12);
null = zeros(1,12);

% Create xls file base for tolerance values:
phioverlapxls(pfilename);

% Calculate tolerances:
for layer = 1:3
    % Read z & r tolerances that make up z-overlap from excel file:
    if layer == 1
       z_tol(1,1:12) = xlsread(zfilename, 'Layer 1', 'A11:L11');
       z_tol(2,1:12) = xlsread(zfilename, 'Layer 1', 'A12:L12');
       r_tol(1,1:12) = xlsread(zfilename, 'Layer 1', 'A14:L14');
       r_tol(2,1:12) = xlsread(zfilename, 'Layer 1', 'A15:L15');
    elseif layer == 2
       z_tol(1,1:12) = xlsread(zfilename, 'Layer 2', 'A11:L11');
       z_tol(2,1:12) = xlsread(zfilename, 'Layer 2', 'A12:L12');
       r_tol(1,1:12) = xlsread(zfilename, 'Layer 2', 'A14:L14');
       r_tol(2,1:12) = xlsread(zfilename, 'Layer 2', 'A15:L15');
    elseif layer == 3
       z_tol(1,1:12) = xlsread(zfilename, 'Layer 3', 'A11:L11');
       z_tol(2,1:12) = xlsread(zfilename, 'Layer 3', 'A12:L12');
       r_tol(1,1:12) = xlsread(zfilename, 'Layer 3', 'A14:L14');
       r_tol(2,1:12) = xlsread(zfilename, 'Layer 3', 'A15:L15');
    end
    
    % Calculate z-overlaps with worst-case module positions:
    for i = 1:12
        z_change_zo(i) = z_tol(1,i);
        r_change_zo(i) = -r_tol(2,i);
        if i > 1
            z_change_zo(i-1) = -z_tol(2,i-1);
            r_change_zo(i-1) = r_tol(1,i-1);
        end 
        z_cov_results_worst = zcov(layer, 'outer', z_change_zo, r_change_zo, tilt_change, 1);
        z_overlap_worst(i) = z_cov_results_worst(12+i);
        z_change_zo = zeros(1,12);
        r_change_zo = zeros(1,12);
    end
    
    % Calculate z-overlaps with nominal module positions:
    z_cov_results_nominal = zcov(layer, 'outer', null, null, null, 1);
    z_overlap_nominal = z_cov_results_nominal(13:24);
    
    if strcmp(phi_results, 'worst') == 1
        tol = pove_tol(layer, z_overlap_worst, z_change, r_change, x_change, x_tilt_change, z_tilt_change, z_ratio, r_ratio, x_ratio, x_tilt_ratio, z_tilt_ratio, limit);
    elseif strcmp(phi_results, 'nominal') == 1
        tol = pove_tol(layer, z_overlap_nominal, z_change, r_change, x_change, x_tilt_change, z_tilt_change, z_ratio, r_ratio, x_ratio, x_tilt_ratio, z_tilt_ratio, limit);
    end
    xlswrite(pfilename, tol(1,1:12)./2, ['Layer ' num2str(layer)], 'B2:M2');
    xlswrite(pfilename, tol(2,1:12)./2, ['Layer ' num2str(layer)], 'B3:M3');
    xlswrite(pfilename, tol(3,1:12)./2, ['Layer ' num2str(layer)], 'B4:M4');
    xlswrite(pfilename, rad2deg(tol(4,1:12)./2), ['Layer ' num2str(layer)], 'B5:M5');
    xlswrite(pfilename, rad2deg(tol(5,1:12)./2), ['Layer ' num2str(layer)], 'B6:M6');
end
%% z-coverage sensitivity analysis
% Quasi-randomize    
modules = 'outer'; % Modules to analyze (tolerances tighter in outer modules across layers, which effectively results in outer module tolerances being
% the actual positional tolerances of the ring)
eval_tol = 'no'; % Use pre-evaluated values for module tolerances?
tol = 5; % If fixed tolerances will be used for analysis, specify magnitude of +- tolerance of z & r movement for modules
if strcmp(eval_tol, 'yes') == 1
    % Read tolerances from tolerance xls file according to layer:
    if layer == 1
       z_tol(1,1:12) = xlsread(filename, 'Layer 1', 'A11:L11');
       z_tol(2,1:12) = xlsread(filename, 'Layer 1', 'A12:L12');
       r_tol(1,1:12) = xlsread(filename, 'Layer 1', 'A14:L14');
       r_tol(2,1:12) = xlsread(filename, 'Layer 1', 'A15:L15');
    elseif layer == 2
       z_tol(1,1:12) = xlsread(filename, 'Layer 2', 'A11:L11');
       z_tol(2,1:12) = xlsread(filename, 'Layer 2', 'A12:L12');
       r_tol(1,1:12) = xlsread(filename, 'Layer 2', 'A14:L14');
       r_tol(2,1:12) = xlsread(filename, 'Layer 2', 'A15:L15');
    elseif layer == 3
       z_tol(1,1:12) = xlsread(filename, 'Layer 3', 'A11:L11');
       z_tol(2,1:12) = xlsread(filename, 'Layer 3', 'A12:L12');
       r_tol(1,1:12) = xlsread(filename, 'Layer 3', 'A14:L14');
       r_tol(2,1:12) = xlsread(filename, 'Layer 3', 'A15:L15');
    end
elseif strcmp (eval_tol, 'no') == 1
    z_tol = repmat(tol,2,12);
    r_tol = repmat(tol,2,12);
end

% Randomize values between tolerance values:
for i=1:12
    z_change(:,i)        =(z_tol(1,i)+z_tol(2,i)).*lhs(runs, 1)-z_tol(2,i);
    r_change(:,i)        =(r_tol(1,i)+r_tol(2,i)).*lhs(runs, 1)-r_tol(2,i);
    tilt_change(:,i)     =(tilt_tol(1,i)+tilt_tol(2,i)).*lhs(runs, 1)-tilt_tol(2,i);
end
% Loop z coverances with randomized module positional:
for i = 1:3
    layer = i;
    for index=1:12
        z_cov_single = zcov_single(layer, modules, index, z_change, r_change, tilt_change, runs);
        Y = z_cov_single;
        % See how change matrix correlate with results vector. Whole
        % change matrix included to account for correlation throughout the
        % tilted section:
        St = plcc([z_change, r_change, tilt_change],Y);
        z_st(index) = sum(St(1:12));
        r_st(index) = sum(St(13:24));
        t_st(index) = sum(St(25:36));
    end
    % Plotting
    figure;
    plot(1:12,z_st); grid on; hold on;
    plot(1:12,r_st); plot(1:12,t_st);
    ylabel('Correlation coefficient'); xlabel('Ring');
    legend('z variation', 'r variation', 'tilt variation');
    if strcmp(eval_tol, 'yes') == 1
        title(['Tolerance sensitivies with pre-evaluated tolerances - Layer ', num2str(layer)]);
    elseif strcmp(eval_tol, 'no') == 1
        title(['z-coverage tolerance sensitivies with fixed tolerances (z & r = +-',num2str(tol), ' mm) - Layer ', num2str(i)]);
    end
end
%% Next two sections can be used for testing:
% Single test run for z-coverage with nominal values
clear
layer = 1;
modules = 'inner';
z_change = zeros(1,12);
r_change = zeros(1,12);
tilt_change = zeros(1,12);
z_cov_result = zcov(layer, modules, z_change, r_change, tilt_change, 1);

%% Single test run for phi-overlap with nominal values
clear
layer = 1;
z_change = zeros(1,12);
r_change = zeros(1,12);
x_change = zeros(1,12);
x_tilt_change = zeros(1,12);
y_tilt_change = zeros(1,12);
z_tilt_change = zeros(1,12);
% Calculate z-overlaps with nominal module positions:
z_change_zo = zeros(1,12);
r_change_zo = zeros(1,12);
tilt_change_zo = zeros(1,12);
z_cov_result_zo = zcov(layer, 'outer', z_change_zo, r_change_zo, tilt_change_zo, 1);
z_overlap = z_cov_result_zo(13:24);
pol = zeros(1,12);

for ring = 1:12
    pol(ring) = phi_overlap(layer, ring, z_change, r_change, x_change, x_tilt_change, y_tilt_change, z_tilt_change, z_overlap, 1);
end