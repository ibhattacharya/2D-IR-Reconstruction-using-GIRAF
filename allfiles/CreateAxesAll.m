%% Create Axes
n = load('ArrayCalibConsts_Shiftfreq_Interval.calib'); % calibration parameters 

% n(1-3); %quadratic fit terms for camera calibration
% n(4); % center wavelength of camera; 
% n(5); % wavelength of up conversion lightn=
% n(6); % freq shift after FT
% n(7); % tau time step in fs

%% Full axes
% w1 - linear based on FFT
w1_min = n(6)+(1/(3*1e-5*1024*n(7))); w1_max = n(6)+(512/(3*1e-5*1024*n(7)));
w1_axis= linspace(w1_min,w1_max,512);

% w3 - nonlinear
w3_axisNL = zeros(1,1024);
for z = 1:1024;
    w3_axisNL(z) = 1*1e7*(1/(n(1)*z^2+n(2)*z+n(3)) - 1/(n(5)));
end

% w3 - linear
w3_axisL = linspace(w3_axisNL(1),w3_axisNL(end),length(w3_axisNL));



