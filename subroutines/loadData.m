%loads data from input file if its in three column format
function [x, y, z] = loadData(dataFile) 
    data = load(dataFile);
    x = data(:, 1);
    y = data(:, 2);
    z = data(:, 3);
end 