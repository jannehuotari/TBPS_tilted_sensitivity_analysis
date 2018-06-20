% This function returns a 5x12 array of tolerances with the following
% format:
% First row = z tolerance
% Second row = r tolerance
% Third row = v tolerance
% Fourth row = beta tolerance
% Fifth row = gamma tolerance

% Layer = layer to analyze (integer)

% z_overlap = z_overlap values to use in analysis (1 x 12 array)

% z, r , x_tilt etc... = 1 x [value of your choosing] arrays that get
% incrementally larger. These increments will be implemented into model of
% phi-overlap to calculate the tolerances.

% z_ratio, r_ratio etc... = ratios that are < 1. When a value is found from
% the incremented arrays where phi-overlap = 0, the value needs to be
% corrected with some ratio to calculate the rest of the tolerances.

% limit = limit for phi-overlap. Functional requirement set at 0, but
% tolerance can be calculated with more tight / lax requirements.

function [tolerances] = pove_tol(layer, z_overlap, z, r, x, x_tilt, z_tilt, z_ratio, r_ratio, x_ratio, x_tilt_ratio, z_tilt_ratio, limit)

array_size = size(z, 2);
% phi_overlap.m was designed foremost for sensitivity analysis, so it takes
% runs x 12 change matrix as input. Change arrays need to be transformed to
% comply:
z_transpose = z';
z_change = repmat(z_transpose, 1, 12);
r_transpose = r';
r_change = repmat(r_transpose, 1, 12);
x_transpose = x';
x_change = repmat(x_transpose, 1, 12);
x_tilt_transpose = x_tilt';
x_tilt_change = repmat(x_tilt_transpose, 1, 12);
z_tilt_transpose = z_tilt';
z_tilt_change = repmat(z_tilt_transpose, 1, 12);

% Initialize a generic zero matrix:
null = zeros(array_size,12);

% Initialize array for tolerances
% First row = z tolerance
% Second row = r tolerance
% Third row = v tolerance
% Fourth row = beta tolerance
% Fifth row = gamma tolerance
tolerances = zeros(5,12);

% Tolerances are determined based on their corresponding dimensions
% influence of phi-overlap. Order: gamma, v, r, beta and z

% Run loops for all rings for all tolerances:

% Gamma shift:
for ring = 1:12
    i = 1;
    while phi_overlap(layer, ring, null, null, null, null, null, z_tilt_change(i,1:12), z_overlap, 1) > limit
        tolerances(5, ring) = z_tilt_change(i,1);
        i = i + 1;
        if i == 100000
            break
        end
    end
end
% Reset i, multiply maximum gamma shift with the given ratio and
% move on to calculate the next tolerance (v), with the calculated gamma shift applied:
tolerances(5, 1:12) = tolerances(5, 1:12) .* z_tilt_ratio;
for ring = 1:12
    i = 1;
    % Figure out at which v shift phi-overlap goes to zero:
    while phi_overlap(layer, ring, null, null, x_change(i,1:12), null, null, tolerances(5, 1:12), z_overlap, 1) > limit
        tolerances(3, ring) = x_change(i,1);
        i = i + 1;
        if i == 100000
            break
        end
    end
end

% r tolerance:
tolerances(3,1:12) = tolerances(3,1:12).* x_ratio;
for ring = 1:12
    i = 1;
    while phi_overlap(layer, ring, null, r_change(i,1:12), tolerances(3,1:12), null, null, tolerances(5,1:12), z_overlap, 1) > limit
        tolerances(2, ring) = r_change(i,1);
        i = i + 1;
        if i == 100000
            break
        end
    end
end

% x_tilt tolerance:
tolerances(2,1:12) = tolerances(2,1:12).* r_ratio;
for ring = 1:12
    i = 1;
    while phi_overlap(layer, ring, null, tolerances(2,1:12), tolerances(3,1:12), x_tilt_change(i,1:12), null, tolerances(5,1:12), z_overlap, 1) > limit
        tolerances(4, ring) = x_tilt_change(i,1);
        i = i + 1;
        if i == 100000
            break
        end
    end
end

% z tolerance (only if layer 1, in 2 and 3 negligible):
tolerances(4,1:12) = tolerances(4,1:12).* x_tilt_ratio;
for ring = 1:12
    i = 1;
    while phi_overlap(layer, ring, z_change(i,1:12), tolerances(2,1:12), tolerances(3,1:12), tolerances(4,1:12), null, tolerances(5,1:12), z_overlap, 1) > limit
        tolerances(1, ring) = z_change(i,1);
        i = i + 1;
        if i == 100000
            break
        end
    end
end
tolerances(1,1:12) = tolerances(1,1:12).*z_ratio;
% y_tilt and z_overlap need not be accounted for, because according to the
% sensitivity analysis, they have a negligible effect on phi-overlap
end
