function [Q,freq_axis,alpha, freq_axis_rad] = calculateQ(path,f,R,fd_check)
% A function that reads different F-D curves from XLSX files and returns
% aveage relaxance.
% [Q,freq_axis] = calculateQ(path,f,R,fd_check)
% Returs relaxance and related frequency axis
% Requires path of the data files, sampling frequency and tip radius and 
% whteher if it is fd curve or stress strain. if FD curve use 'fd' if not
% it will assume stress strain curve.




% First we need to import the data. Data import can be done according to
% data files we have. In this case we have forces and indentations in excel
% files.



%Here I am going to take a look at fd_check. If its fd, add contact
%mechancis.
fd_flag = strcmp(fd_check,'fd');
power =  fd_flag  * 3/2 + (~fd_flag);

%% This includes data import part. You should arrange this according to 
%your own data structure

%Getting the file names in the given path
fileList = dir([path '*.xlsx']);

%Some initialization
alpha = 0;

for i = 1:length(fileList)
    
    %Reading the data
    file_name = fileList(i).name;
    tmp = xlsread([path file_name]);
    
    

    %input data
    data(i).indent = tmp(:,1)' ;%- tmp(1,1)';
    data(i).force = tmp(:,2)' ;%- tmp(1,2)';

    %some moving average
    mov_win = ceil(length(tmp(:,1)) * 0.08);
    mov_win = mov_win + (mov_win == 0);
    
    
    data(i).indent = movmean(data(i).indent, mov_win);
    data(i).force = movmean(data(i).force, mov_win);

    
    %finding optipzed alpha for the modified Fourier transfom
    tmp_oa_f = optimizedAlpha(data(i).force);
    tmp_oa_h = optimizedAlpha(data(i).indent.^(power));
    %This (3/2) is coming from contact mechanics
    
    alpha = max([tmp_oa_h tmp_oa_f alpha]);
    lengths(i) = length(data(i).indent);


end

alpha = alpha * 10;
%In here we are going to average the spectrums to get rid of noise. Since
%every spectrum have slightly different number of elements, it is not
%stright forward. We are going to use interpolation. We are going to
%interpolate to the shortest force - distance 
% This is goiung to be the freuqnecy axis for it in the unit circle
freq_axis_rad = linspace(-pi, pi, min(lengths));





%%

fourier_length = min(lengths);
spectrum_h = zeros(length(fileList) ,fourier_length);
spectrum_f = zeros(length(fileList) ,fourier_length);
for i = 1:length(fileList)
    %Calculating the modified spectrum for force and indentation^(3/2) 
    %as the contact mechanics dictate 10.31181/rme200102156m
    spectrum_f(i,:) = fftshift( fft(data(i).force ...
        .* exp(-alpha .* (1:(length(data(i).force)))),...
        fourier_length)) * 2 / f;
    
    spectrum_h(i,:) = fftshift( fft(data(i).indent.^(power) ...
        .* exp(-alpha .* (1:(length(data(i).indent)))), ...
        fourier_length) ) * 2 / f;
    

end


if fd_flag
    Q = (sum(spectrum_f) ./ sum(spectrum_h)) * 9 / (16*sqrt(R));
else
    Q = (sum(spectrum_f) ./ sum(spectrum_h));
end


%here we calculate the relaxance and retardance.
freq_axis = freq_axis_rad * (f/2) / pi;





end
