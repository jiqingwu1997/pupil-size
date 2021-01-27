function VideoON = TT_DetectVideoON(Filebase, nChs, Fs)

% 10 Oct 2018
% this function will determine the onset of pupil monitoring
% and make .video.evt file WITHOUT making any sync file

% ASSUMPTION
% A single sync pulse appears within first 30 sec in Ch#1.
% the pulse must be the first

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% params
if nargin<3; Fs = 20000; end
Seg = 30; % loading data segment in sec
SyncCh = 1;
maxVolts = 5;
smpsPerVolt = double(intmax('int16')/maxVolts);
HamWin = 8;
LookAhead = 4*(Fs/1000);
syncThresh = 0.15;         % exclude events on sync chan (trigCh) less than 0.1V
syncRefrac = 4*(Fs/1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load first 
tmp =  LoadDatSeg([Filebase, '.dat'], nChs, 1, Fs*Seg);
sync = tmp(SyncCh,:)/smpsPerVolt;
clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% event detection
b = diff(hamming(HamWin),1); b = b-mean(b);   % diff of smoothed sig
dEeg = Filter0(b,sync);                          % find event times & amplitudes
minima = LocalMinima(-abs(dEeg), syncRefrac, -syncThresh)-1; % To find LocalMinima, see /ken's code/General/

evtTimes = minima((minima>LookAhead & minima<length(dEeg)-LookAhead));
evtAmps = sync(evtTimes+LookAhead) - sync(evtTimes-LookAhead);

MasterEvtAmps = evtAmps;
MasterEvtTimes = evtTimes/(Fs/1000);  % ms-based

%% visualization
Time=(1:length(sync))/(Fs/1000);
figure(1);plot(Time/1000, sync,'k-');hold on;
plot(MasterEvtTimes/1000, MasterEvtAmps, 'k.');hold off;
xlabel('time (sec)'); ylabel('amp (V)');
axis([0 Seg -1 1.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% video on
VideoON = MasterEvtTimes(1); % in ms
% create file
dlmwrite([Filebase,'.video.evt'],[VideoON 1],'-append', 'delimiter', '\t', 'precision', '%.2f');

% update fig
figure(1);hold on;
plot([VideoON VideoON]/1000,[-5 5],'r:');hold off;

