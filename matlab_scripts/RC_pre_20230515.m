% MEP Recruitment Curve Script
% Last Edit: May-15-2023
clear

addpath('/rri_disks/artemis/meltzer_lab/jed/toolboxes/eeglab2020_0')

%%%%%%%%%%%%%%%%%%%%%%%%%%% Script Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output Parameters (1=on, 0=off)
plot_trials = 1; % matlab plot
return_trials = 1;
return_avg_trials = 1;

% Trial rejection params
rejection = 1; % EEGLAB trial rejection (1=on, 0=off)
trials_to_reject = []; % Trials to automatically skip (default=[])

% Set window for identifying MEP max and min values. In line 8, 
% set plot_frames = 1 to scroll through each trial and determine the
% optimal start_look (usually 0.01 [10ms] - 0.03 [30ms]) and end_look 
% (usually around 0.08 [80ms]) times
mep_start_time = 0.02;
mep_end_time = 0.04;

sample_rate = 5000;

%%%%%%%%%%%%%%%%%%%%%%%%%% Read in Raw Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ask user to select file
[file_name, path_name] = uigetfile('*.csv', 'Select the RC data');
path = fullfile(path_name, file_name); % file path
protocol = 'mepRC'; % Protocol (a part of the file name unique to protocol)

subID = strsplit(path_name, "SUBJID_"); % pulls out subject ID from pathname 
subID = strsplit(subID{2}, "/");
subID = subID{1};

% Check that correct file was selected, end program if not
if not(contains(path, protocol))
    fprintf('File does not match protocol: %s \n', protocol)
    return
end

% Read in raw data
raw_data = readmatrix(path, 'NumHeaderLines', 8);
raw_data = raw_data(:, 1:end-1); % remove last column to prevent import errors
clipped_data = rmmissing(raw_data(2:end, :)); % clip time column

% Define hemisphere specific parameters for loading channel data
if contains(path, 'Rhem')
    channel = 1;
    hem = 'Rhem';
elseif contains(path, 'Lhem')
    channel = 2;
    hem = 'Lhem';
end

% Define useful variables relating to the dataset
frames = raw_data(:, 1); % Frame column data
states = raw_data(:, 2); % State column data
chans = raw_data(2:end, 3); % Channel # column data
pulse_times = raw_data(:, 6); % Pulse time column

% Get unique values (remove duplicates)
unique_frames = rmmissing(unique(frames)); 
unique_states = rmmissing(unique(states));
unique_chans = rmmissing(unique(chans));

% Get number of unique values for each column
num_frames = length(unique_frames);
num_states = length(unique_states);
num_chans = length(unique_chans);

% Length of the channel data (from start time to end time)
num_epochs = length(clipped_data(:, 10:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filter Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array to store re-formatted channel data
chan_data = zeros(num_chans, num_epochs, num_frames);

% Fill chan_data with eeg data, contains all channels and frames
for i = 1:length(unique_chans)
    chan = unique_chans(i);
    chan_inds = find(chans == chan);
    chan_data(chan, :, :) = raw_data(chan_inds+1, 10:end)';
end

% Import channel data into EEGLAB
eeglab nogui
EEG = pop_importdata('dataformat', 'array', 'nbchan', 0, 'data', 'chan_data', 'srate', 5000, 'pnts', 0, 'xmin', 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname','emg', 'gui', 'off'); 
EEG = eeg_checkset(EEG);

% Apply notch filter
% set 'revfilt' to 1, inverts filter (from bandpass to notch filter)
% {default=0 (bandpass)}
EEG = pop_eegfiltnew(EEG, 'locutoff', 58.5, 'hicutoff', 61.5,  'filtorder', 16500, 'revfilt', 1, 'plotfreqz', 0);

% Apply bandpass filter
EEG = pop_eegfiltnew(EEG, 'locutoff', 25.5, 'hicutoff', 999.5, 'filtorder', 16500, 'plotfreqz', 0);

% Store eeglab data
cleaned_data = EEG.data;

% New matrix to store re-formatted channel data
clean_formatted = zeros(size(cleaned_data, 3) * num_chans, length(raw_data));
% Add channel, trial, pulse information
clean_formatted(:, 1:9) = raw_data(2:end, 1:9);

% Indexes for each channel (skip by num chans)
formatting_inds = 1:num_chans:size(cleaned_data, 3)*num_chans;

% Add channels from EEGLAB dataset to re-fromatted
for i = 1:size(EEG.data, 3)
    clean_formatted(num_chans*i-1:num_chans*i, 10:end) = EEG.data(:, :, i);
end

% Store reformatted datasets
clipped_data = clean_formatted;
raw_data = vertcat(raw_data(1, :), clean_formatted);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Find MEPs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create variables to store MEP data
trial_data = cell(0, 0);

% Run through each frame and find MEP measurements
for i = 1:length(unique_frames)
    frame = unique_frames(i);
    frame_inds = find(frames == frame); % Row indexes for current frame #
    State = states(frame_inds); State = State(1); % Current state
    
    % Load channel data (ch1=EMG, ch2=pulse)
    if channel == 1
        ch1 = raw_data(frame_inds(1), 10:end);
        ch2 = raw_data(frame_inds(2), 10:end);
    elseif channel == 2
        ch1 = raw_data(frame_inds(2), 10:end);
        ch2 = raw_data(frame_inds(1), 10:end);
    end
    
    time = raw_data(1, 10:end); % Get time data
    pulse_time = raw_data(frame_inds(channel), 6); % TMS pulse time
    pulse_length = raw_data(frame_inds(channel), 7); % Pulse duration
    pulse_sample = sample_rate * pulse_time; % TMS pulse time (# samples)

    start_look = round(sample_rate * (pulse_time + mep_start_time));
    end_look = round(sample_rate * (pulse_time + mep_end_time));

    % Fist find the max and min values from 10-ms post pulse to end of a reasonable window
    [maxVAL, maxIND] = max(ch1(start_look:end_look));                    
    maxIND = maxIND + start_look;                                                            

    % find the minimum MEP value            
    [minVAL, minIND] = min(ch1(start_look:end_look));
    minIND = start_look + minIND;
    
    % calculate the peak-to-trough amplitude
    mepAMP = maxVAL - minVAL;

    % If the minimum peak comes before the maximum, take latency as the
    % time of the minimum. Otherwise (time of max is before time of min
    % peak), take latency as the time of the maximum peak (usual scenario)
    if minIND < maxIND
        mep_latency = 1000 * (minIND / sample_rate - pulse_time);
    else
        mep_latency = 1000 * (maxIND / sample_rate - pulse_time);
    end
    
    frame_vals = {frame, State, pulse_time, mepAMP, maxVAL, minVAL,...
        maxIND * 1000 / sample_rate, minIND * 1000 / sample_rate,...
        mep_latency};
    
    % Calculate size of frame vals
    if i == 1
        trial_data = cell(0, length(frame_vals));
    end
    
    % Store trial data
    trial_data(frame, :) = frame_vals;
    
    %%% Plot MEP %%%
    if plot_trials
        % Plot EMG data
        plot(time(625:1125) * 1000, ch1(625:1125) , '-b');
        title(['RC-', hem, ' State: ', num2str(State), ', Frame: ', num2str(frame)])
        
        % Axis labels
        xlabel('Time (ms)')
        ylabel('Filtered EMG (mV)')
        
        % Plot pulse + mep
        xline(pulse_time * 1000, '-b', {'Pulse'})
        xline(time(maxIND) * 1000,'-g', {'MEP Max'})
        xline(time(minIND) * 1000,'-r', {'MEP Min'})
        
        xline(start_look * 1000 / sample_rate, '-b', {'Start look'})
        xline(end_look * 1000 / sample_rate, '-b', {'End look'})
        
        pause
    end
end

rejected_trials = []; % to store rejected trial indexes

if and(rejection, isempty(trials_to_reject))
    % Array to store re-formatted channel data
    chan_data = zeros(num_chans, num_epochs, num_frames);

    % Fill chan_data with eeg data, contains all channels and frames
    for i = 1:length(unique_chans)
        chan = unique_chans(i);
        chan_inds = find(chans == chan);
        chan_data(chan, :, :) = raw_data(chan_inds+1, 10:end)';
    end
    
    mep_pulse_times = cell2mat(trial_data(:, 3)) / 1000;
    mep_max = cell2mat(trial_data(:, 7)) / 1000;
    mep_min = cell2mat(trial_data(:, 8)) / 1000;
    
    event_info = [rmmissing(unique(frames)), mep_pulse_times / sample_rate, mep_max, mep_min];
    event_info = num2cell(event_info);
    
    % Import 'EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_importdata('dataformat', 'array', 'nbchan', 0, 'data', 'chan_data', 'srate', 5000, 'pnts', 0, 'xmin', 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname','emg', 'gui', 'off'); 
    EEG = eeg_checkset(EEG);
    
    % Add event markers (MEP/SP start/end, pulse times)
    EEG = pop_importepoch(EEG, event_info, {'epoch', 'pulse', 'mep_max', 'mep_min'}, 'latencyfields', {'pulse', 'mep_max', 'mep_min'}, 'timeunit', 1, 'headerlines', 0);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );
    
    % Plot data in EEGLAB
    eeglab redraw
    raw_plot = pop_eegplot(EEG, 1, 1);
    
    all_EEG = EEG.data;
    
    % Button that pauses script until user presses button
    message = sprintf('Press OK once EEGLAB is closed.');
    uiwait(msgbox(message)); 
    
    % Extract cleaned data from EEGLAB (will not contain rejected trials)
    cleaned_data = EEG.data;
    
    % New dataframe to store re-formatted data
    clean_formatted = zeros(size(cleaned_data, 3) * num_chans, length(raw_data));
    
    % Indexes for for each frame (skip by num_chans to avoid indexing same
    % frame
    all_inds = 1:num_chans:size(cleaned_data, 3)*num_chans;
    clipped_inds = [];
    rejected_trials = [];
    
    % Find rejected trial indexes
    for i = 1:size(all_EEG, 3)
        found = 0;
        for n = 1:size(cleaned_data, 3)
            if all_EEG(:, 1:30, i) == cleaned_data(:, 1:30, n)
                found = n;
            end
        end
        
        if not(found)
            rejected_trials(end+1) = i;
        end
    end
    
    % Reformat EEGLAB data to script format
    for i = 1:size(cleaned_data, 3) % For each frame
        frame_vals = cleaned_data(1, 1:30, i); % Take first 30 vals of data
        start_frame_ind = 0;

        % Find corresponding frame_vals in raw_data to get frame # using
        % the first 30 values to compare
        for row = 2:num_chans:size(raw_data, 1)
            if frame_vals == single(raw_data(row, 10:39))
                start_frame_ind = row;
                break
            end
        end
        
        if start_frame_ind > 0
            % Load frame data and store it
            frame_inds = start_frame_ind:start_frame_ind+num_chans-1;
            frame_data = raw_data(frame_inds, :);
            clean_inds = all_inds(i):all_inds(i)+num_chans-1;
            clean_formatted(clean_inds, :) = frame_data;
            
            % Store index and update to next
            clipped_inds(end+1) = start_frame_ind;
        end
    end
    
    clipped_inds = clipped_inds / num_chans;
    
    % Re-define MEP data with rejected trials removed
    trial_data = trial_data(clipped_inds, :);
    
elseif not(isempty(trials_to_reject))
    clipped_inds = unique_frames;
    rejected_trials = trials_to_reject;
    
    for i = trials_to_reject
        clipped_inds = clipped_inds(find(clipped_inds~=i));
    end
    
    % Re-define MEP data with rejected trials removed
    trial_data = trial_data(clipped_inds, :);
end

% Return average data to .csv file
if return_avg_trials
    % Cell to store average values
    avg_vals = cell(0, 10);
    
    for i = 1:num_states
        % Indexes for current state
        state_inds = find(cell2mat(trial_data(:, 2)) == i);
        pulse = trial_data(state_inds(1), 3);
        mepAMP = mean(nonzeros(cell2mat(trial_data(state_inds, 4))));
        maxVAL = mean(nonzeros(cell2mat(trial_data(state_inds, 5))));
        minVAL = mean(nonzeros(cell2mat(trial_data(state_inds, 6))));
        maxT = mean(nonzeros(cell2mat(trial_data(state_inds, 7))));
        minT = mean(nonzeros(cell2mat(trial_data(state_inds, 8))));
        mepLAT = mean(nonzeros(cell2mat(trial_data(state_inds, 9))));
        
        % Calculate number of rejected trials
        n_rejected = length(find(states == i)) / num_chans - length(state_inds);
        n_trials = length(state_inds);
        
        avg_val = {i, pulse, mepAMP, maxVAL, minVAL, maxT, minT, mepLAT,...
            n_trials, n_rejected};
        avg_vals(end+1, :) = avg_val;
    end
    
    list_of_rejects = cell(num_states, 1);
    for state_ = 1:num_states
        list_of_rejects{state_} = "";
    end
    
    % Add list of rejected trial inds
    if rejected_trials
        for trial = rejected_trials
            state_ = states(trial * num_chans);
            list_of_rejects{state_} = append(list_of_rejects{state_}, " " + num2str(trial));
        end
    end
    avg_vals(:, end+1) = list_of_rejects;
    
    % Output file name
    file_name = [subID, '_pre-', protocol, '-', hem, '-avg.csv'];
    
    % Write output file
    output_file = array2table(avg_vals, 'VariableNames',...
        {'State', 'Pulse Time (ms)', 'MEP AMP (mV)', 'MEP Max (mV)',...
        'MEP Min (mV)', 'MEP Max T. (ms)', 'MEP Min T. (ms)',...
        'MEP Latency (ms)', '# Trials', '# Rejected', 'Rejected Trials'});
    writetable(output_file, file_name)
end

% If specified to return individual frame data, return csv file
if return_trials    
    file_name = [subID, '_pre-', protocol, '-', hem, '-trials.csv']; % File name
    
    output_file = array2table(trial_data, 'VariableNames',...
        {'Frame', 'State', 'Pulse', 'MEP AMP (mV)', 'MEP Max (mV)',...
        'MEP Min (mV)', 'MEP Max T. (ms)', 'MEP Min T. (ms)',...
        'Latency (ms)'});
    writetable(output_file, file_name) % Output file
end
