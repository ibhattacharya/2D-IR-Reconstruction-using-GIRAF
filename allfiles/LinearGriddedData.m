%% Linear Gridded 2D Data, Humston 2017
% The w3 axis is slightly nonlinear. This code is to be used to interpolate the data 
% onto a linearized grid. 
% The inputs are data and the outputs of CreateAxes code. The inputs are saved with 
% new names to designate nonlinearity (dataNL and w3_axisNL) and the outputs are saved
% with the original names (data and w3_axis), this way this piece of code can be inserted
% into any other code without disrupting it, and also without losing the origianl matrices.

%% Run with data and calib files already loaded and CreateAxes run
s= size(data);

%The original data is saved as dataNL 
dataNL = data;
clear data;

% Generate linear w3 axis - determine the query points for sampling over
% the range of w3
%3_axisL = linspace(w3_axisNL(w3Trunc_low),w3_axisNL(w3Trunc_high),s(1));

% Data to be interpolated
data = zeros(size(dataNL)); % Holds linearized gridded data
for i=1:size(dataNL,2);
    
    v = dataNL(:,i);
    % Interpolate the function at the query points
    data(:,i) = interp1(w3_axisNL,v,w3_axisL,'spline');
end
