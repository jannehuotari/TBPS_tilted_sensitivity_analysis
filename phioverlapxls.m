% This function produces an xls file that serves as a base to fill in
% tolerance values related to phi-overlap

function phioverlapxls(filename)
xlswrite(filename, cellstr('Layer 1'), 'Layer 1', 'A1:A1');
xlswrite(filename, cellstr('z-tolerance'), 'Layer 1', 'A2:A2');
xlswrite(filename, cellstr('r-tolerance'), 'Layer 1', 'A3:A3');
xlswrite(filename, cellstr('v-tolerance'), 'Layer 1', 'A4:A4');
xlswrite(filename, cellstr('beta-tolerance'), 'Layer 1', 'A5:A5');
xlswrite(filename, cellstr('gamma-tolerance'), 'Layer 1', 'A6:A6');

xlswrite(filename, cellstr('Layer 2'), 'Layer 2', 'A1:A1');
xlswrite(filename, cellstr('z-tolerance'), 'Layer 2', 'A2:A2');
xlswrite(filename, cellstr('r-tolerance'), 'Layer 2', 'A3:A3');
xlswrite(filename, cellstr('v-tolerance'), 'Layer 2', 'A4:A4');
xlswrite(filename, cellstr('beta-tolerance'), 'Layer 2', 'A5:A5');
xlswrite(filename, cellstr('gamma-tolerance'), 'Layer 2', 'A6:A6');

xlswrite(filename, cellstr('Layer 3'), 'Layer 3', 'A1:A1');
xlswrite(filename, cellstr('z-tolerance'), 'Layer 3', 'A2:A2');
xlswrite(filename, cellstr('r-tolerance'), 'Layer 3', 'A3:A3');
xlswrite(filename, cellstr('v-tolerance'), 'Layer 3', 'A4:A4');
xlswrite(filename, cellstr('beta-tolerance'), 'Layer 3', 'A5:A5');
xlswrite(filename, cellstr('gamma-tolerance'), 'Layer 3', 'A6:A6');
end