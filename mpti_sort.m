function [r_same, r_diff] = mpti_sort(tab)

%% 1: Conversion of filename into condition ID (same/different and SOA)
% L/R = presentation side for visual stimulus
% S/D = same / different spatial condition. Where the audio appeared?
% Number =   0    1    2    3    4    5   6  7  8  9   10  11  12  13  14
%          -350 -300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300 350
%          temporal condition. minus means audio leading conditions
%
% e.g. from LS7.mp4 to Left presentation of flash | same spatial condition | 8th temporal condition (0 to 14)

cond = tab.movieFile;
resp = tab.rating_keykeys;
cond(1:12) = []; % deletes first 12 rows (the practice phase)
resp(1:12) = [];
mat = zeros(length(cond),3); % final matrix to build

for j=1:length(cond)
    % spatial congruency - same / diff. Side is not considered
    if startsWith(cond{j},'LS') || startsWith(cond{j},'RS') 
        mat(j,1) = 1; % same position
    else
        mat(j,1) = 2; % diff position
    end
    
    % temporal condition (SOA)
    % Converting from string to double, it tries on the 3rd and 4th chars
    % of the string. They could be a double-digit number (else) or a single-digit with the dot for the extension
    % included in the conversion. In this case, it shows a NaN.
    if isnan(str2double(cond{j}([3 4])))         % if the string is not a number (e.g. 3_),
        mat(j,2) = 1 + str2double(cond{j}(3));   % only convert the first digit
    else
        mat(j,2) = 1 + str2double(cond{j}([3 4]));
    end
    
    % response
    % Modified on 28/07/2021: distinction between certainty levels
    % no more grouping
    mat(j,3) = resp(j);
    
end

%% 2: Each response is added to an array for that Spatial and Temp condition
% heaps for same position condition
t1 = []; t2 = []; t3 = []; t4 = []; t5 = []; t6 = []; t7 = []; 
t8 = []; t9 = []; t10 = []; t11 = []; t12 = []; t13 = []; t14 = []; t15 = [];
% heaps for different position condition (14 = 1, ...)
t16 = []; t17 = []; t18 = []; t19 = []; t20 = []; t21 = [];
t22 = []; t23 = []; t24 = []; t25 = []; t26 = []; t27 = []; t28 = []; t29 = []; t30 = [];
% (There is surely a better way)

for i = 1:length(mat)
    if mat(i,1) == 1 % same position
        switch mat(i,2) % based on SOA, add the current response to the corresponding slot
            case 1, t1 = horzcat(t1,mat(i,3)); 
            case 2, t2 = horzcat(t2,mat(i,3)); 
            case 3, t3 = horzcat(t3,mat(i,3)); 
            case 4, t4 = horzcat(t4,mat(i,3));
            case 5, t5 = horzcat(t5,mat(i,3));
            case 6, t6 = horzcat(t6,mat(i,3));
            case 7, t7 = horzcat(t7,mat(i,3));
            case 8, t8 = horzcat(t8,mat(i,3));
            case 9, t9 = horzcat(t9,mat(i,3));
            case 10, t10 = horzcat(t10,mat(i,3));
            case 11, t11 = horzcat(t11,mat(i,3));
            case 12, t12 = horzcat(t12,mat(i,3));
            case 13, t13 = horzcat(t13,mat(i,3));
            case 14, t14 = horzcat(t14,mat(i,3));
            case 15, t15 = horzcat(t15,mat(i,3));
        end
    else
        switch mat(i,2)
            case 1, t16 = horzcat(t16,mat(i,3));
            case 2, t17 = horzcat(t17,mat(i,3));
            case 3, t18 = horzcat(t18,mat(i,3));
            case 4, t19 = horzcat(t19,mat(i,3));
            case 5, t20 = horzcat(t20,mat(i,3));
            case 6, t21 = horzcat(t21,mat(i,3));
            case 7, t22 = horzcat(t22,mat(i,3));
            case 8, t23 = horzcat(t23,mat(i,3));
            case 9, t24 = horzcat(t24,mat(i,3));
            case 10, t25 = horzcat(t25,mat(i,3));
            case 11, t26 = horzcat(t26,mat(i,3));
            case 12, t27 = horzcat(t27,mat(i,3));
            case 13, t28 = horzcat(t28,mat(i,3));
            case 14, t29 = horzcat(t29,mat(i,3));
            case 15, t30 = horzcat(t30,mat(i,3));
        end
    end
end

%% 3: Data cleaning
%
% Matrix format (every coloumn is a SOA): 
%   - Row1: # 5 answers (surely simultaneous)
%   - Row2: # 4 answers (maybe simultaneous)
%   - Row3: # 3 answers (non lo so)
%   - Row4: # 2 answers (maybe non simultaneous)
%   - Row5: # 1 answers (surely non simultaneous)
%   - Row6: Number of responses (different from 3)
%   - Row7: #5 / (#5 + #1) aka "Sure" rate
%   - Row8: #4 / (#4 + #2) aka "Maybe" rate

r_same = zeros(8,15); 
r_diff = zeros(8,15); 

for k = 1:15
    % same position
    eval(['r_same(1,k) = sum(t' num2str(k) '(:) == 5);']); 
    eval(['r_same(2,k) = sum(t' num2str(k) '(:) == 4);']);
    eval(['r_same(3,k) = sum(t' num2str(k) '(:) == 3);']);
    eval(['r_same(4,k) = sum(t' num2str(k) '(:) == 2);']);
    eval(['r_same(5,k) = sum(t' num2str(k) '(:) == 1);']);
    r_same(6,k) = r_same(1,k) + r_same(2,k) + r_same(4,k) + r_same(5,k);
    r_same(7,k) = r_same(1,k) / (r_same(1,k) + r_same(5,k));
    r_same(8,k) = r_same(2,k) / (r_same(2,k) + r_same(4,k)); 
    
    % diff position
    eval(['r_diff(1,k) = sum(t' num2str(k) '(:) == 5);']); 
    eval(['r_diff(2,k) = sum(t' num2str(k) '(:) == 4);']);
    eval(['r_diff(3,k) = sum(t' num2str(k) '(:) == 3);']);
    eval(['r_diff(4,k) = sum(t' num2str(k) '(:) == 2);']);
    eval(['r_diff(5,k) = sum(t' num2str(k) '(:) == 1);']);
    r_diff(6,k) = r_diff(1,k) + r_diff(2,k) + r_diff(4,k) + r_diff(5,k);
    r_diff(7,k) = r_diff(1,k) / (r_diff(1,k) + r_diff(5,k));
    r_diff(8,k) = r_diff(2,k) / (r_diff(2,k) + r_diff(4,k));
end

end