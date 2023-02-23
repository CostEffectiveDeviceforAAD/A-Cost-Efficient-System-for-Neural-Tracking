%% Before Experiment - Window Size 15 sec. / for Real

clear all

%%%%%%%%%%%%%%%%%%%%%%%%% Set-up Directory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set directory
cd 'C:\Users\Seung-Cheol Baek\AppData\Roaming\Neuroscan\Curry 8\Matlab'

% display the working directory
disp(pwd)


%%%%%%%%%%%%%%%%%%%%%%%%% ??? LSL Part ??? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% +) import LSL library


%%%%%%%%%%%%%%%%%%%%%% Set-up Filtering Params %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ====================== Original Data Info ============================

% original sampling rate of eeg
orig_sr = 1000; % in Hz

% ===================== Filtering Params ==============================

% nyquist frequency
nyquist = orig_sr/2;

% lower & high filter bound
lower_filter_bound = 2; % Hz
upper_filter_bound = 8; % Hz

% transition width for filter construction
transition_width   = 0.25;

revfilt = 0; % low & band pass filter, 1 for high-pass

% ============ Create filter weights as in pop_eegnewfilt ================

edgeArray = sort([lower_filter_bound upper_filter_bound]);

% Max stop-band width
maxTBWArray = edgeArray; % Band-/highpass
if revfilt == 0 % Band-/lowpass
    maxTBWArray(end) = nyquist - edgeArray(end);
end
maxDf = min(maxTBWArray);

% Transition band width and filter order
df = min([max([edgeArray(1) * transition_width 2]) maxDf]);

% filter_order
filter_order = 3.3 / (lower_filter_bound / orig_sr); % Hamming window
filter_order = ceil(filter_order / 2) * 2; % Filter order must be even.

% Passband edge to cutoff (transition band center; -6 dB)
dfArray = {df, [-df, df]; -df, [df, -df]};
cutoffArray = edgeArray + dfArray{revfilt + 1, length(edgeArray)} / 2;

% Window
winArray = windows('hamming', filter_order + 1);

% Filter coefficients
global filterweights
filterweights = firws(filter_order, cutoffArray / nyquist, winArray);


%%%%%%%%%%%%%%%%%%%%%%%%% Figure Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************** You should check up values assigned to "k" (3 or 8)

figure(1)
% Plot for BOTH Decoder
subplot(2,1,1)
global h1; global h2;
h1=plot(0, 0, 'LineWidth', 1.3);
hold on
h2=plot(0, 0, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*8],'ylim',[-0.5 0.5])
for k = 1:8
    if k ~= 1
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend([h1, h2], {'Audio From Left', 'Audio From Right'})
title(' Both Decoder ', 'FontWeight', 'bold')

% Plot for Biased Decoder
subplot(2,1,2)
global h3; global h4;
h3=plot(0, 0, 'LineWidth', 1.3);
hold on
h4=plot(0, 0, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*8],'ylim',[-0.5 0.5])
for k = 1:8
    if k ~= 1
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend([h3, h4], {'Audio From Left', 'Audio From Right'})
title(' Biased Decoder ', 'FontWeight', 'bold')



%%%%%%%%%%%%%%%%%%%%%%% Delete Local Vars. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

