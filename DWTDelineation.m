% This file converts ECG data into numerical arras

clear;
close all;
clc;
% Include directory of your files.
files = dir("C:\Users\fanga\OneDrive\Documents\Fall 2023\elec 594\THC Docs\TCH Data\Dataset 1\JET data\**\*h5");
filesList = []; % This Variable will store the names of all files as a list.
for i=1:length(files)
    filesList = [filesList ; strcat(files(i).folder,'\',files(i).name)];
end
%filesList(1,114:131) name index for SINUS data
%filesList(1,112:129) name index for JET data
%filesList(12,:) = [];
[k, m] = size(filesList);
nameList = [];
sizeList = [];
% Perform iteration to extract recordings from each file. 
for i=1:k
    filename = filesList(i,:);
    
    info = h5info(filename);
    h5disp(info.Filename);


    HR = h5read(info.Filename,[info.Name 'PARM_HR']);
    ECG1 = h5read(info.Filename,[info.Name 'GE_WAVE_ECG_1_ID']);
    ECG2 = h5read(info.Filename,[info.Name 'GE_WAVE_ECG_2_ID']);
    ECG2 = reshape(ECG2, [1,length(ECG2)]);
    time = h5read(info.Filename,[info.Name 'time']);

    %dwt1 = computeDWT(ECG2);
    signal = ECG2;
    fs = 1/mean(diff(time));
    % bandpass filter 0.5-50
    filtSignal = bandpass(signal, [0.5, 50], fs);
    %butter third order 5hz high pass.
    [b,a] = butter(3, 5/(fs/2), 'high');
    filtSignal = filter(b, a, filtSignal);

    array = computeDWT(signal); % Compute DWT
    [modLines1, max2, min2] = findModLinesQRS(array); % Process onset and offset of R wave

    [newLines, groups] = processBeats(modLines1, HR, 240); % make sure peaks are properly spaced.


    limits = delineator(newLines, array(4,:), fs);
    
    heartbeat = zeros(length(limits) - 1, 100);
    
    for j = 1 : length(limits) - 1
        beg = limits(j);
        fin = limits(j+1);
        record = filtSignal(beg:fin-1);
        samPoints = linspace(1,100, length(record));
        x = linspace(1,100,100);
        heartbeat(j,:) = interp1(samPoints, record, x);
    end

    [rows,cols] = size(heartbeat);
    
    %Normalize heartbeats as probabilities.
    for n=1:rows
        heartbeat(n,:) = heartbeat(n,:) - min(heartbeat(n,:));
        heartbeat(n,:) = heartbeat(n,:)/ sum(heartbeat(n,:));
    end
    writematrix(heartbeat, strcat(filename(112:129),'.csv')); 
    
    nameList = [nameList;filename(112:129)];
    sizeList = [sizeList; size(heartbeat)];
    

end
% writematrix(nameList, 'nameList.csv');
% writematrix(sizeList, 'sizeList.csv');

%Define functions
function DWT = computeDWT(signal)
%    This function will take in a raw ECG signal and will return the 4 scales of the DWT 
%    using the algorithm described in Martinez et, al.
%    The output will be a 2D matrix of size 4xN with N as the length of original ECG signal
    
%     Note: the filters used in this function assume the signal is sampled at 250 Hz. Most signals from
%         this dataset are sampled at 240Hz. Therefore, resampling the filters might be a possible improvement
%     

    % Create the FIR of low pass and high pass filter.
    
    x = linspace(-10,9,20);
    
    h = (1/8)*(deltaDiscrete(x+2) + 3*deltaDiscrete(x+1) + 3*deltaDiscrete(x) + deltaDiscrete(x-1) );
    g = 2*( deltaDiscrete(x+1) - deltaDiscrete(x) );
    
    N = length(signal);
    scale1 = conv(signal, g);
    scale1 = scale1(find(scale1~=0,1): find(scale1~=0,1) + N - 1);
    temp = conv(signal, h);
    scale2 = conv(temp,g);
    scale2 = scale2(find(scale2~=0,1): find(scale2~=0,1) + N - 1);
    temp = conv(temp, h);
    scale3 = conv(temp,g);
    scale3 = scale3(find(scale3~=0,1): find(scale3~=0,1) + N - 1);
    temp = conv(temp,h);
    scale4 = conv(temp,g);
    scale4 = scale4(find(scale4~=0,1): find(scale4~=0,1) + N - 1);
    

    DWT = [scale1; scale2; scale3; scale4];
    
end



function [modLines, maxima2, minima2] = findModLinesQRS(array)
% This function will scan an array of DWT and will return the modulus
% lines.
    N = length(array);
    
    maxima2 = zeros(1,N);
    minima2 = zeros(1,N);
    W = 65536;
    % Find Thresholds   
    k = 1;
    l = 1;
    for i = 1:W:N
        epsi4 = rms(array(4,i:min(N, i + W - 1)));
        epsi3 = rms(array(3,i:min(N, i + W - 1)));
        epsi2 = rms(array(2,i:min(N, i + W - 1)));
        %epsi1 = rms(array(1,i:min(N, i + W - 1)));

        [~, maxima] = findpeaks(array(4,i:min(N, i + W - 1)), 'MinPeakHeight', epsi4);
        [~, minima] = findpeaks(-array(4,i:min(N, i + W - 1)), 'MinPeakHeight', epsi4);
        
        maxima = matchPeaks(maxima, array(3,i:min(N, i + W - 1)), epsi3);
        maxima = matchPeaks(maxima, array(2,i:min(N, i + W - 1)), epsi2) + i;
        %maxima = matchPeaks(maxima, array(1,i:min(N, i + W - 1)), epsi1) + i;
        maxima2(k:k + length(maxima) - 1) = maxima;
        k = k + length(maxima);
        
        minima = matchPeaks(minima, -array(3,i:min(N, i + W - 1)), epsi3);
        minima = matchPeaks(minima, -array(2,i:min(N, i + W - 1)), epsi2) + i;
        %minima = matchPeaks(minima, -array(1,i:min(N, i + W - 1)), epsi1) + i;
        minima2(l:l + length(minima) - 1) = minima;
        l = l + length(minima);
        
    end 
    maxima2 = maxima2(1:k-1);
    minima2 = minima2(1:l-1);
    
    
    
    modLines = findModLines(maxima2, minima2);
end

function newPeaks = matchPeaks(idx, array, epsi)
    
    newPeaks = zeros(1,length(idx));
    j = 1;
    for i = 1:length(idx)
        k = idx(i);
        if k <=5
            continue
        end
        if k+5 > length(array)
            continue
        end
        logical = array(k-5:k+5) > epsi;
        if any(logical == true)
            [~, maxIdx] = max(array(k-5:k+5));
            offset = maxIdx - 6;
            newPeaks(j) =  k + offset;
            j = j + 1;
        end
    end
    newPeaks = newPeaks(1:j-1);
        

end

function modLines = findModLines(maxima, minima)
    modLines = zeros(2, length(maxima));
    k = 1;
    j = 1;
    for i=1:length(maxima)
        idx = maxima(i);
        mini = inf;
        pair = 0;
        if j > length(minima)
            break
        end
        
        while minima(j) < idx || abs( minima(j) - idx ) <=10
            
            dist = abs(minima(j) - idx);
            if dist <= 10 && dist <= mini
                mini = dist;
                pair = minima(j);                
            end
            j = j + 1;
            if j > length(minima)
                break
            end
        end
        if pair ~= 0
            modLines(1,k) = min(pair, idx);
            modLines(2,k) = max(pair, idx);
            k = k + 1;
        end
    end
    modLines = modLines(:,1:k-1);
    
end

function [newMod, groups] = processBeats(modulusArray, heartRate, fs)
    
    newMod = modulusArray;
    
    sortHR = sort(heartRate, "descend");
    if length(sortHR) >= 5
        maxHR = sortHR(5);
    else
        maxHR = max(heartRate);
    end

    epsi = 60*fs/maxHR - 15;
    
    t = 1;
    while t < 10
    
        diffArray = diff(newMod(1,:));

        ids = find(diffArray < epsi);


        groups = zeros(2, length(ids));

        i = 1;
        k = 1;
        while i < length(ids)
            groups(1,k) = ids(i);
            while i < length(ids) && ids(i+1) == ids(i) + 1
                groups(2,k) = ids(i+1);
                i = i + 1;
            end
            i = i + 1;
            k = k + 1;
        end
        
        groups = groups(:, 1:k-1);
        for i = 1:size(groups,2)
            lastPeak = groups(1,i);
            idx = groups(1,i) + 1;
            lastIdx = groups(2,i);
             if lastIdx == 0
                 newMod(1,idx) = 0;

             else
                while idx <= lastIdx
                    dist = modulusArray(1,idx) - modulusArray(1,lastPeak);
                    if dist < epsi
                        newMod(1,idx) = 0;
                    else
                        lastPeak = idx;

                    end
                    idx = idx + 1;
                end
            end
        end
        newIdx = newMod(1,:) ~= 0;
        newMod = newMod(:,newIdx); 
        if isempty ( find ( diff(newMod(1,:)) < epsi, 1 ) ) == 1
            break
        end
        t = t + 1;
    
    end
    
    
end

function limits = delineator(modLines, scale4, fs)
    
    % Generate array of thresholds.
%     W = 65536;
%     N = length(scale4);
%     thresholds = zeros(1,ceil(N/W));
%     for i=1:W:N
%         thresholds(i) = 0.15*rms(scale4(i: min(i + W - 1, N)));
%     end
    
    limits = zeros(1,length(modLines));
    k = 1;
    offSet = floor((0.04*fs));
    
    Tduration = ceil(0.1 * fs);
    
    for i = 1:length(modLines)-1
        
        
        
        if (modLines(1, i + 1) - modLines(2,i) > 0.7 * fs)
            windowT = floor(0.55 * fs);
        else
            windowT = floor( 0.65 * (modLines(1, i + 1) - modLines(1,i) ));
        end
        
        
        
        endi = modLines(2,i);
        window = scale4(endi + offSet: endi + offSet + windowT);        
        epsi = rms(window);
        
        
        maxima = find(window > epsi);
        minima = find(-window > epsi);
        
        if isempty(maxima) == 0 && isempty(minima) == 0
            min1 = -1;
            endT = inf;
            for t=1:length(minima)
                m = minima(t);
                l = 1;
                while l <= length(maxima)
                    dist = abs(m - maxima(l));
                    
                    if dist < Tduration && dist > min1
                        endT = max(m, maxima(l));
                        min1 = dist;
                    end
                    l = l + 1;
                end
            end  
            
            if endT == inf
                endT = max(find( window == max(window) ), find( window == min(window) ));
            end
            
            zeroCrossing = find ( abs(diff(sign(window(endT:length(window)) ) ) ) > 0, 1, 'first');
            if isempty(zeroCrossing) == 0
                endT = endT + zeroCrossing;
            end
            limits(1,k) = endT + endi + offSet - 2;
            
            k = k + 1;
            
        else
            limits(1,k) = endi + windowT;
            
            k = k + 1;
        end
        
    end

limits = limits(1:k-1);
end

function delta = deltaDiscrete(x)
% This function will simulate a discrete delta dirac function where the output of an array will be zero except 
% at point x=0 where the output will be 1.
    delta = zeros(1,length(x));
    for i=1:length(x)
        if x(i) == 0
            delta(i) = 1;
        end
    end
    
end


