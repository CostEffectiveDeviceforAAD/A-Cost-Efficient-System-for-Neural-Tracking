%% Neuroscan Simulation
%
% written by Seung-Cheol Baek
% Update Date: Mon. 13th July, 2020
% 


%% Figure Preparation

f_h1=figure(1);
clf

f_h2=figure(2);
clf
set(f_h2, 'Color', 'w');

global h1
h1=plot(0, 0, 'LineWidth', 1.3);
hold on
global h2
h2=plot(0, 0, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*6],'ylim',[-0.5 0.5])
for k = 1:6
    if k ~= 1
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend([h1, h2], {'Audio From Left', 'Audio From Right'})


f_h3=figure(3);
clf
set(f_h3, 'Color', 'w');

global h3
h3=plot(0, 0, 'LineWidth', 1.3);
hold on
global h4
h4=plot(0, 0, 'LineWidth', 1.3);
grid on
set(gca,'xlim',[1 46*6],'ylim',[-0.5 0.5])
for k = 1:6
    if k ~= 1
        plot( [(k-1)*46+1 (k-1)*46+1], get(gca,'ylim'), 'k--', 'LineWidth', 0.5)
    end
end
xlabel('Time (sec.)'), ylabel('Correlation Coefficients')
legend([h3, h4], {'Audio From Left', 'Audio From Right'})


%% LSL Preparation

% instantiate the library
disp('Loading the library...');
global lib
lib = lsl_loadlib();

% resolve a stream...
disp('Resolving an EEG stream...');
eeg = {};
while isempty(eeg)
    eeg = lsl_resolve_byprop(lib,'type','EEG');  end

disp('Resolving an Marker stream...');
marker = {};
while isempty(marker)
    marker = lsl_resolve_byprop(lib,'type','Marker'); end

% create a new inlet
disp('Opening an inlet...');
global eeg_inlet
eeg_inlet = lsl_inlet(eeg{1});
global marker_inlet
marker_inlet = lsl_inlet(marker{1});


%% Simulation


% initialize i
i = 0;

while true
   
   tic;
    
   % update count
   i = i + 1;
   
   % File embedded to Neuroscan
   ACQ1OnlineAAD; 
   
   % Break
   if i > 690
       break
   end
   
   t=toc;
   
   pause(1-t+0.0163)
   
end
