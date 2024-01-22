% Short Interval Cortical Inhibition Script
% Last Edit: Mar-16-2023
clear 


%%%%%%%%%%%%%%%%%%%%%%%%% Output Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output parameters (1=On, 0=off)
plot_frames = 1; % Plot individual frames (matlab plot)
return_trials = 1; % Return individual trial data in csv
return_states = 1; % Return averaged trial data based on state

% Trial rejection params
rejection = 1; % EEGLAB trial rejection (1=on, 0=off)
trials_to_reject = []; % Trials to automatically skip (default=[])

% Set window for identifying MEP max and min values. In line 8, 
% set plot_frames = 1 to scroll through each trial and determine the
% optimal start_look (usually 0.01 [10ms] - 0.03 [30ms]) and end_look 
% (usually around 0.08 [80ms]) times
mep_start_time = 0.02;
mep_end_time = 0.045;

sample_rate = 5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read in Raw Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ask user to select file
[file_name, path_name] = uigetfile('*.csv', 'Select the SICI data');
path = fullfile(path_name, file_name);
protocol = 'sici';

subID = strsplit(path_name, "SUBJID_"); % pulls out subject ID from pathname 
subID = strsplit(subID{2}, "/");
subID = subID{1};

% Define channel/hem. of interest (contralateral hand?)
if contains(path, 'Rhem')
    channel = 1;
    hem = 'Rhem';
elseif contains(path, 'Lhem')
    channel = 2;
    hem = 'Lhem';
end

% Check that correct file was selected, end program if not
if not(contains(path, protocol))
    fprintf('File does not match protocol: %s \n', protocol)
    return
end

% Read in data
raw_data = readmatrix(path, 'NumHeaderLines', 8);
raw_data = raw_data(:, 1:end-1);
clipped_data = rmmissing(raw_data(2:end, :));

% Define useful variables relating to the dataset
frames = rmmissing(raw_data(:, 1)); % Frame column data
states = rmmissing(raw_data(:, 2)); % State column data
chans = rmmissing(raw_data(2:end, 3)); % Channel # column data
pulse_times = rmmissing(raw_data(:, 6)); % Pulse time column
time = rmmissing(raw_data(1, 10:end));

% Get unique values (remove duplicates)
unique_frames = unique(frames); 
unique_states = unique(states);
unique_chans = unique(chans);

% Get number of unique values for each column
num_frames = length(unique_frames);
num_states = length(unique_states);
num_chans = length(unique_chans);

% Length of the channel data (from start time to end time)
num_samples = length(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filter Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Array to store re-formatted channel data
chan_data = zeros(num_chans, num_samples, num_frames);

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

% Matrix to store re-formatted data
clean_formatted = zeros(size(cleaned_data, 3) * num_chans, length(raw_data));

% Add channel, trial, pulse information
clean_formatted(:, 1:9) = raw_data(2:end, 1:9);

% Indexes for each channel
formatting_inds = 1:num_chans:size(cleaned_data, 3)*num_chans;

% Add channels from EEGLAB dataset to re-fromatted
for i = 1:size(EEG.data, 3)
    clean_formatted(num_chans*i-1:num_chans*i, 10:end) = EEG.data(:, :, i);
end

clipped_data = clean_formatted;
raw_data = vertcat(raw_data(1, :), clean_formatted);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Find SICI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0 placeholder, rename variable on first iteration to calculate size of
% data to store
trial_data = cell(0, 0);

for i = unique_frames'
    % Define frame specific variables
    frame = i;
    frame_inds = find(frames == frame);
    frame_data = clipped_data(frame_inds, :);
    
    % Load in channel data
    ch1 = frame_data(1, 1:10);
    ch2 = frame_data(2, 1:10);
    
    state = ch2(2);
    tag = ch2(4);
    
    % Assign hemisphere specific variables 
    if channel == 1
        data = frame_data(1, 10:end);
        num_pulses = ch1(5);
        p1_time = ch1(6);
        p1_length = ch1(7);
        p2_time = ch1(8);
        p2_length = ch1(9);
        hemisphere = 'Right';
    elseif channel == 2
        data = frame_data(2, 10:end);
        num_pulses = ch2(5);
        p1_time = ch2(6);
        p1_length = ch2(7);
        p2_time = ch2(8);
        p2_length = ch2(9);
        hemisphere = 'Left';
    end
    
    % determine the conditioning and test pulses based on the pulse timing
    % values (i.e., it is a test only if the second pulse time is 0)
    if p2_time > 0
        CT_pulse_times = [p1_time, p2_time];
    else
        CT_pulse_times = [p1_time, 0];
    end
    
    % MEP Channel data
    rect_data = abs(data);
    diff_data = horzcat(0, diff(data));
    rect_diff = abs(diff_data);
      
    % Variables to store MEP values
    CTmaxVAL = [0,0];
    CTmaxIND = [0,0];
    CTminVAL = [0,0];
    CTminIND = [0,0];
    CTmepAMP = [0,0];
    CTmepLAT = [0,0];
    
    % Only look for T pulse (2nd) MEP as C pulse MEP doesn't occur for SICI
    if state == 1
        pulse_ind = 1;
    else
        pulse_ind = 2;
    end
    
    % Set window for identifying MEP max and min values. In line 8, 
    % set plot_frames = 1 to scroll through each trial and determine the
    % optimal start_look (usually 0.01 [10ms] - 0.03 [30ms]) and end_look 
    % (usually around 0.08 [80ms]) times
    test_pulse = CT_pulse_times(pulse_ind);
    start_look = round((test_pulse + mep_start_time) * sample_rate);
    end_look = round((test_pulse + mep_end_time) * sample_rate);

    % Find max and min values, as well as time of max and min peaks (ind), 
    % between 10ms post-stimulus and the specified end look time
    [v_max, ind_max] = max(data(start_look:end_look));
    CTmaxVAL(pulse_ind) = v_max;
    CTmaxIND(pulse_ind) = ind_max + start_look;

    [v_min, ind_min] = min(data(start_look:end_look));
    CTminVAL(pulse_ind) = v_min;
    CTminIND(pulse_ind) = ind_min + start_look;
    
    % If the minimum peak comes before the maximum, take latency as the
    % time of the minimum. Otherwise (time of max is before time of min
    % peak), take latency as the time of the maximum peak (usual scenario)
    if ind_min < ind_max
        CTmepLAT(pulse_ind) = 1000 * (CTminIND(pulse_ind) / sample_rate - test_pulse); 
    else
        CTmepLAT(pulse_ind) = 1000 * (CTmaxIND(pulse_ind) / sample_rate - test_pulse);
    end

    % MEP Amplitude   
    CTmepAMP(pulse_ind) = CTmaxVAL(pulse_ind) - CTminVAL(pulse_ind);

    %%%%%%%%%%%%%%%%%%%%%%%%% Store Trial Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Variable for trial data
    if state == 1
        frame_vals = {frame, state, CT_pulse_times(2)*1000, CT_pulse_times(1)*1000,...
            CTmepAMP(1), CTmaxVAL(1), CTminVAL(1), CTmaxIND(1) * 1000 / sample_rate,... 
            CTminIND(1) * 1000 / sample_rate,... 
            CTmepLAT(1)};
    else
        frame_vals = {frame, state, CT_pulse_times(1)*1000, CT_pulse_times(2)*1000,...
            CTmepAMP(2), CTmaxVAL(2), CTminVAL(2), CTmaxIND(2) * 1000 / sample_rate,...
            CTminIND(2) * 1000 / sample_rate,...
            CTmepLAT(2)};
    end

    % Calculate size of frame vals
    if i == 1
        trial_data = cell(0, length(frame_vals));
    end
    
    % Store trial data
    trial_data(frame, :) = frame_vals;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% Plot Trial Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if plot_frames
        % Plot EMG data
        plot(time * 1000, data, 'b'); hold on
        plot(time * 1000, data, 'b');
        title(['SICI-', hem, 'Frame: ', num2str(frame), 'State: ', num2str(state)])
        
        % Axis labels
        xlabel('Time (ms)')
        if channel == 1
            ylabel('Left Hand (mV)')
        elseif channel == 2
            ylabel('Right Hand (mV)')
        end
        
        % Plot pulse times and meps
        if state == 1
            xline(CT_pulse_times(1) * 1000, '-c', {'Test'})
            xline(CTminIND(1) * 1000 / sample_rate, 'r', {'Min'})
            xline(CTmaxIND(1) * 1000 / sample_rate, 'g', {'Max'})

        else
            xline(CT_pulse_times(1) * 1000, '-c', {'Condition'})
            xline(CT_pulse_times(2) * 1000, 'c', {'Test'})
            xline(CTminIND(2) * 1000 / sample_rate, 'r', {'Min'})
            xline(CTmaxIND(2) * 1000 / sample_rate, 'g', {'Max'})
        end
        
        xline(start_look * 1000 / sample_rate, '-b', {'Start look'})
        xline(end_look * 1000 / sample_rate, '-b', {'End look'})
        
        hold off
        pause
    end
end

rejected_trials = [];

% EEGLAB Trial rejection and plotting
if and(rejection, isempty(trials_to_reject))
    % Array to store re-formatted channel data
    chan_data = zeros(num_chans, num_samples, num_frames);

    % Fill chan_data with eeg data, contains all channels and frames
    for i = 1:length(unique_chans)
        chan = unique_chans(i);
        chan_inds = find(chans == chan);
        chan_data(chan, :, :) = raw_data(chan_inds+1, 10:end)';
    end
    
    % Create cell to store event info
    pulse_times1 = cell2mat(trial_data(:, 3)) / 1000;
    pulse_times2 = cell2mat(trial_data(:, 4)) / 1000;
    mep_max = cell2mat(trial_data(:, 8)) / 1000;
    mep_min = cell2mat(trial_data(:, 9)) / 1000;
    
    event_info = [rmmissing(unique(frames)), pulse_times1, pulse_times2,...
        mep_max, mep_min];
    event_info = num2cell(event_info);

    % Import channel data into EEGLAB
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_importdata('dataformat', 'array', 'nbchan', 0, 'data', 'chan_data', 'srate', 5000, 'pnts', 0, 'xmin', 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname','emg', 'gui', 'off'); 
    EEG = eeg_checkset(EEG);
    
    % Add event markers (MEP max/min, pulse times)
    EEG = pop_importepoch(EEG, event_info, {'epoch', 'pulse1', 'pulse2', 'mep_max', 'mep_min'}, 'latencyfields', {'pulse1', 'pulse2', 'mep_max', 'mep_min'}, 'timeunit', 1, 'headerlines', 0);
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
    
    % Remove rejected trials from raw_data, cycling each frame to check if
    % cleaned data is in raw data (all frames)
    for i = 1:size(cleaned_data, 3)
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
    
    % Re-define SICI data with rejected trials removed
    old_trial_data = trial_data;
    trial_data = trial_data(clipped_inds, :);
    
elseif not(isempty(trials_to_reject))
    clipped_inds = unique_frames;
    rejected_trials = trials_to_reject;
    
    for i = trials_to_reject
        clipped_inds = clipped_inds(find(clipped_inds~=i));
    end
    
    old_trial_data = trial_data;
    trial_data = trial_data(clipped_inds, :);
end
    
% Return individual trials as csv
if return_trials
    file_name = [subID, '_pre', protocol, '-', hem, '-trials.csv']; % File name
   
    output_file = array2table(trial_data, 'VariableNames',...
        {'Frame', 'State', 'C. Pulse Time (ms)', 'T. Pulse Time (ms)',...
        'MEP AMP. (mV)', 'MEP Max (mV)', 'MEP Min (mV)',...
        'MEP Max T. (ms)', 'MEP Min T. (ms)', 'MEP Latency T. (ms)'});
    writetable(output_file, file_name) % Output file
end

% Return averaged state
if return_states
    file_name = [subID, '_pre', protocol, '-', hem, '-avg.csv'];
    avg_vals = cell(0, 11);
    
    for State = 1:num_states
        state_inds = find(cell2mat(trial_data(:, 2)) == State);
        pulse1 = trial_data(state_inds(1), 3);
        pulse2 = trial_data(state_inds(1), 4);
        mepAMPs = mean(nonzeros(cell2mat(trial_data(state_inds, 5))));
        maxVALs = mean(nonzeros(cell2mat(trial_data(state_inds, 6))));
        minVALs = mean(cell2mat(trial_data(state_inds, 7)));
        minINDs = mean(cell2mat(trial_data(state_inds, 8)));
        maxINDs = mean(cell2mat(trial_data(state_inds, 9)));
        mepLATs = mean(cell2mat(trial_data(state_inds, 10)));
        
        n_trials = length(find(cell2mat(trial_data(:, 2)) == State));
        
        if rejection
            all_trials = length(find(cell2mat(old_trial_data(:, 2)) == State));
            n_rejected = abs(all_trials - length(state_inds));
        else
            n_rejected = 0;
        end
        
        avg_vals(end+1, :) = {State, pulse1, pulse2, mepAMPs,...
            maxVALs, minVALs, minINDs, maxINDs, mepLATs, n_trials, n_rejected};
    end
    
    list_of_rejects = cell(num_states, 1);
    for state_ = 1:num_states
        list_of_rejects{state_} = "";
    end
    
    if rejected_trials
        for trial = rejected_trials
            state_ = states(trial * num_chans);
            list_of_rejects{state_} = append(list_of_rejects{state_}, " " + num2str(trial));
        end
    end
    
    avg_vals(:, end+1) = list_of_rejects;
    
    output_file = cell2table(avg_vals, 'VariableNames', {'State',...
        'C. Pulse Time (ms)', 'T. Pulse Time (ms)', 'MEP AMP. (mV)',...
        'MEP Max (mV)', 'MEP Min (mV)', 'MEP Max T. (ms)',...
        'MEP Min T. (ms)', 'MEP Latency T. (ms)', '# Trials', '# Rejected', 'Rejected Trials'});
    writetable(output_file, file_name);
end
