% All Rights Reserved.
% Copyright (c) 2025 Md Mayenul Islam.
% Rajshahi University of Engineering & Technology, Bangladesh.
% Email: mayenulislam15@eee.ruet.ac.bd.

% Citation for this work:
% [1] M. M. Islam, M. Abdul Motin, S. Kabir, and D. Kumar, “Poincaré image 
% analysis of short-term electrocardiogram for  detecting atrial fibrillation,” 
% in IEEE Journal of Biomedical and Health Informatics, accepted, 2025, 
% doi: 10.1109/JBHI.2025.3562778R.

% The entire process will require the following MATLAB Toolboxes:
            %1. wfdb toolbox
            %2. ECG-kit toolbox
            %2. CurveLab toolbox
            %3. Empirical Wavelet Transforms(https://matlab.mathworks.com/open/fileexchange/v1?id=42141)
% Label = 0 for Normal, 1 for AF, 3 for PACs/PVCS, Later label 3 is replaced by 0 considering PACs/PVCs as normal
% It's a demo code using a record from CPSC TR-I. It will provide hints to implement the paper

clc;
clear all;
SubID = 'data_3_1';
DataPath = 'CPSC21 Dataset\';
[signal, Fs, tm] = rdsamp([DataPath, SubID]);
[SigInfo] = wfdbdesc([DataPath, SubID]);
ECGSignal = signal(:, strcmp({SigInfo.Description}, 'II'));
[ann, ann_type, ann_subtype, chan, num, comments] = rdann([DataPath, SubID], 'atr');
%ECGkit toolbox will be needed for finding P locations 
load("p_locations.mat");
r_locations = ann;
% Replace 'None' with '(AFIB': Persistent AF are denoted by 'None
for entry = 1:length(comments)
    if strcmp(comments{entry}, 'None')
        comments{entry} = '(AFIB';
    end
end

unique_chars = unique(ann_type); 
valid_comments = comments(cellfun(@ischar, comments));
unique_strings = unique(valid_comments);
num_beats = length(ann);
labels = cell(num_beats, 1);

last_label = '(N';
for i_beats = 1:num_beats
    if ~isempty(comments{i_beats})
        last_label = comments{i_beats};
    end
    labels{i_beats} = last_label;
end

% Segmentation
si = 0; % segment number
fs = Fs;
second = 10;
val = ECGSignal;
valval = length(val)/(second*Fs);
prev_si = 0;
nonaf = 0;
af = 0;
pac = 0;
pvc = 0;
count=0;
nan =0;

for si = 1:fix(valval)

    l1 = second*fs*(si-1); % segment lower index
    l2 = second*fs*((si-1)+1); % segment upper index
    in_segment = ann >= l1+1 & ann <= l2;
    ann_in_segment = ann(in_segment);
    % Check if the segment has fewer than 5 R-peak positions
    if isempty(ann_in_segment) || length(ann_in_segment) < 5
        continue;  
    end
    segment_labels = labels(in_segment); 
    segment_charaters = ann_type(in_segment);
    unique_labels = unique(segment_labels); % checking the types of labels in the segment
    unique_charaters = unique(segment_charaters);
    if any(strcmp(unique_labels,'(N')) && any(strcmp(unique_labels,'(AFIB'))
        if sum(strcmp(segment_labels, '(AFIB'))>= 0.6*length(segment_labels)
            new_label = 1; %consider as AF
            af = af+1;
        else
            new_label = NaN;
        end
        
    elseif any(strcmp(unique_labels,'(N')) && any(strcmp(unique_labels,'(SBR'))
        new_label = NaN;%3;
        %pac=pac+1;
    elseif any(strcmp(unique_labels,'(N')) || any(strcmp(unique_labels,'(AB')) || any(strcmp(unique_labels,'(B')) || any(strcmp(unique_labels,'(VT')) || any(strcmp(unique_labels,'(SVTA')) %||any(strcmp(unique_labels,'(SBR'))  %strcmp(unique_labels{1}, 'N')
        % Check corresponding ann_type characters for 'A' or 'V'
        if sum(strcmp(segment_labels, '(N'))>= 0.8*length(segment_labels)
            if any(ismember(segment_charaters,['A', 'V'])) %sum(ismember(segment_charaters, ['A', 'V'])) >= 0.1*length(segment_labels)
                new_label = 3;  % Set new_label to 3 if 'A' or 'V' exists
                pac = pac+1;
            elseif any(ismember(segment_charaters,'Q'))
                new_label = NaN;
            elseif sum(ismember(segment_charaters, 'N'))==length(segment_labels)
                new_label = 0;  % Set new_label to 0 if neither exists
                nonaf = nonaf+1;
        end
        else
            new_label = NaN;
        end
    elseif strcmp(unique_labels, '(AFIB')%strcmp(unique_labels{:}, 'AFIB')
        if sum(strcmp(segment_labels, '(AFIB'))>= 0.6*length(segment_labels)
            new_label = 1; %consider as AF
            af = af+1;
        else
            new_label = NaN;
        end
    else
        new_label = NaN;  % Handle unknown labels, if any
    end
    %% RR & dRR determination for the corresponding segment
    if new_label==0 || new_label==1 || new_label==3
        
        r_peak = ann(in_segment)./Fs;
        RR = diff(r_peak);
        dRR = diff(RR);
        RR_Interval = dRR;

        ppl=0;
        ppij=ppi;
        ppil=1;
        [ppif]=0;
    
        for ppil=1:length(ppij)
            if ppij(ppil)>l1+1 && ppij(ppil)<=l2
                ppl=ppl+1;
                ppif(ppl)=ppij(ppil);
            end
        end
        p_peak_length = length(ppif);



        if new_label==0 || new_label==1 || new_label==3
            poincare_image_generation;
            class = new_label;
            if class ==3 %Convert PAC/PVC labels as normal
                class = 0;
            end
            r_p_ratio = p_peak_length/length(r_peak);

            % check if meet the P/R or dRR threshold
            if r_p_ratio>0.7 
                indicator = 1;
            elseif 0.3< r_p_ratio <0.7
                lower_threshold = 0.1;
                upper_threshold = 0.7;
                thp = 0.4;
                a = numel(RR_Interval);
                pir = sum(abs(RR_Interval) >= lower_threshold & abs(RR_Interval) <= upper_threshold) / a;
                if thp > pir && sum(abs(RR_Interval) > upper_threshold) / a > 0
                    indicator = 1; 
                elseif sum(abs(RR_Interval) < lower_threshold) / a >= 0.9
                    indicator = 1; 
                else
                    indicator = 0;
                end
            else
                indicator = 0;
            end



            feature_vector(si,:) = [EA1 EA2 EH1 EH2 EV1 EV2 ED1 ED2 Entr_A1 Entr_A2 Entr_H1 Entr_H2 Entr_V1 Entr_V2 ...
              Entr_D1 Entr_D2 Nor_Entr_A1 Nor_Entr_A2 Nor_Entr_H1 Nor_Entr_H2 Nor_Entr_V1 Nor_Entr_V2 ...
              Nor_Entr_D1 Nor_Entr_D2 E_ewt1 E_ewtf Entr_ewt1 Entr_ewtf Entr_ewt1_nor Entr_ewtf_nor Energy_Appcoeff ...
              Entropy_Appcoeff Entropy_nor_Appcoeff Energy_original Entropy_original Entropy_nor_original Energy_detailcoeff ...
              Entropy_detailcoeff Entropy_nor_detailcoeff  r_p_ratio indicator class];
        else
            continue
        end

   else
        continue
   end


end