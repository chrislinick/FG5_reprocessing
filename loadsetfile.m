function setTable = loadsetfile(filename, dataLines)
%LOADSET Imports all data from a g9 set file into a MATLAB table
%
% INPUT
% path_to_file = the name of file or path to file
%
% OUTPUT
% setTable = all columns of the g9 set file (times and statistics for each set)
% Note: A column of MATLAB datetimes of each set is appended to the table
% (first column), to make plotting and data comparison in MATLAB easier.
%
% EXAMPLE:
% setTable = loadset('S:/Documents/my_set_file.txt');

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [5, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 27);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Set", "Time", "DOY", "Year", "Gravity", "Sigma", "Error", "Uncert", "Tide", "Load", "Baro", "Polar", "Transfer", "Refxo", "Tilt", "Diffraction", "SelfAttract", "Temp", "Pres", "Chan5", "Chan6", "Chan7", "Chan8", "Chan9", "Chan10", "Accept", "Reject"];
opts.VariableTypes = ["double", "string", "string", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Time", "DOY", "Year"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Time", "DOY", "Year"], "EmptyFieldRule", "auto");

% Import the data
setTable = readtable(filename, opts);

% Convert the times to MATLAB datetime format and add to table as column
t_str = setTable.Year + ' ' + setTable.DOY + ' ' + setTable.Time;
t = datetime(t_str,'InputFormat','uuu DDD HH:mm:ss');

setTable = addvars(setTable,t,'before','Set');

end