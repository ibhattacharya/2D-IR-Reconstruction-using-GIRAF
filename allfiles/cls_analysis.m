function [cls] = cls_analysis(data_final,opts)

zoom = data_final(opts.w3_low:opts.w3_high,opts.w1_low:opts.w1_high); %w3 axis comes first because it is rows then columns. The range inputs comes from CreateAxes.m

[M,I] = min(zoom); %Taking minimum for each w1 slice (minimum value from each column matrix)
[m,i] = min(M); %i is a scalar, the location (along w1) of the center of the peak

centpeak = i + opts.w1_low;
centpeakcm = opts.n(6)+((centpeak)/(3*1e-5*1024*opts.n(7))); %Check position of the 'peak' to make sure it is the desired

for l = i - 2: i + 2
    temp = zoom(:,l); %A slice along w1
    temp = movingmean(temp,15);% Taking a moving mean to smooth the daClusterta
    [peakloc,peakmag] = peakfinder(temp,(max(temp)-min(temp))/2,0,-1,true,true); %finds location (index) in w3 of peak for the lth w1 slice
    [y,Y] = min(peakmag); truepeak = peakloc(Y); %Finds the largest peak
    slope(l-i+3,1) = l; slope(l-i+3,2) = truepeak; %Saves cls data, column1 is w1 (in data pt) and column 2 is w3 peak location (in data pt)
end

% Convert point number to wavenumber
w1 = slope(:,1) + opts.w1_low; %w1 is column 1, and adding in the w1 offset, in order to convert from data pt to cm-1
w1 = opts.n(6)+((w1)/(3*1e-5*1024*opts.n(7))); %Converts w1 values from data pt to cm-1.

w3 = slope(:,2) + opts.w3_low + opts.w3Trunc_low; %w3 is column 2, and adding in w3 offset in order to convert from data pt to cm-1
w3_linefit = polyfit(1:length(opts.w3_axisL),opts.w3_axisL,1); %convert w3 linear axis to a function
w3 = w3*w3_linefit(1) + w3_linefit(2);
%Fit array is now two columns, w1 and corresponding peak location along w3.

clsfit = polyfit(w1,w3,1); %Fit to a line
cls = clsfit(1);
