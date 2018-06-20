% This function produces an xls file that serves as a base to fill in
% tolerance values related to z-coverage

function zcovxls(filename)
xlswrite(filename, cellstr('Layer 1'), 'Layer 1', 'A1:A1');
xlswrite(filename, cellstr('Inner'), 'Layer 1', 'A2:A2');
xlswrite(filename, cellstr('z-tolerance'), 'Layer 1', 'A3:A3');
xlswrite(filename, cellstr('r-tolerance'), 'Layer 1', 'A6:A6');
xlswrite(filename, cellstr('Outer'), 'Layer 1', 'A9:A9');
xlswrite(filename, cellstr('z-tolerance'), 'Layer 1', 'A10:A10');
xlswrite(filename, cellstr('r-tolerance'), 'Layer 1', 'A13:A13');

xlswrite(filename, cellstr('Layer 2'), 'Layer 2', 'A1:A1');
xlswrite(filename, cellstr('Inner'), 'Layer 2', 'A2:A2');
xlswrite(filename, cellstr('z-tolerance'), 'Layer 2', 'A3:A3');
xlswrite(filename, cellstr('r-tolerance'), 'Layer 2', 'A6:A6');
xlswrite(filename, cellstr('Outer'), 'Layer 2', 'A9:A9');
xlswrite(filename, cellstr('z-tolerance'), 'Layer 2', 'A10:A10');
xlswrite(filename, cellstr('r-tolerance'), 'Layer 2', 'A13:A13');

xlswrite(filename, cellstr('Layer 3'), 'Layer 3', 'A1:A1');
xlswrite(filename, cellstr('Inner'), 'Layer 3', 'A2:A2');
xlswrite(filename, cellstr('z-tolerance'), 'Layer 3', 'A3:A3');
xlswrite(filename, cellstr('r-tolerance'), 'Layer 3', 'A6:A6');
xlswrite(filename, cellstr('Outer'), 'Layer 3', 'A9:A9');
xlswrite(filename, cellstr('z-tolerance'), 'Layer 3', 'A10:A10');
xlswrite(filename, cellstr('r-tolerance'), 'Layer 3', 'A13:A13');
end