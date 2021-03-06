function output = deapimagcov(filename, dataLines)
%IMPORTFILE Import data from a text file
%  DEAPIMAGCONV = IMPORTFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as a table.
%
%  DEAPIMAGCONV = IMPORTFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  deapimagconv = importfile("/scratch/user/griffith/OoMA-omniscient/src/data/deap_imag_conv.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 09-Sep-2020 15:44:16

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Real", "Imaginary", "RealDEAP", "IMAGDEAP"];
opts.VariableTypes = ["string", "string", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Real", "Imaginary"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Real", "Imaginary"], "EmptyFieldRule", "auto");

% Import the data
deapimagconv = readtable(filename, opts);

end