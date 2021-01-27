function TT_PupilStateAnalysis(Filebase)

% 10 Oct 2018
% this function correlates between pupil diameter and sleep scoring

% requirements
% 1. .pupil.mat (OUT.mat, .L, .R, .D, .t; params; H5list)
% 2. manual sleep scoring\sleepscoring_results.mat (Score)
% 3. .video.evt (on in ms, id)

close all;

% parameter 
sBin = 4; % bin size for sleep scoring (in sec)
Win = 30; % look-bin window (in sec)

%% loading files
load(filebase['C:\Users\lenovo\Downloads\1111\ProcessedData\exp83']; filename=[filebase,'exp83.pupil.mat']);
load('manual sleep scoring\sleepscoring_results.mat');
tmp = load(['C:\Users\lenovo\Downloads\1111\ProcessedData\exp83', 'exp83.video.evt']);

%% extracing distance information
Distance = OUT.D;
Distance = (Distance - mean(Distance))/std(Distance); % z-scored
%% filtering
params.cutF = 0.5;
Distance = lowpass(Distance,params.cutF,params.fps);

%% 
L = OUT.L; R = OUT.R;
Mid = [R(:,1) - L(:,1), abs(R(:,2) - L(:,2))]; % eye center
mMid = mean(Mid); % mean center
dMid = sqrt((Mid(:,1)-mMid(1)).^2+(Mid(:,2)-mMid(2)).^2); % deviation 
dMid = (dMid - mean(dMid))/std(dMid); % z-scored

fps = params.fps;
vON = tmp(1,1);

%% update .pupil.mat file
pT = [1:length(Distance)]/fps + vON/1000; % time vector of pupil data
save(['C:\Users\lenovo\Downloads\1111\ProcessedData\exp83', 'exp83.pupil.mat'], 'pT', '-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sDur = length(Score)*sBin/60; % in min
pDur = (length(Distance)/fps+vON/1000)/60; % in min
dDur = (sDur - pDur)*60; % in sec
fprintf(['Score:', num2str(sDur), 'min\n']);
fprintf(['Distance:', num2str(pDur), 'min\n']);
fprintf(['Score - Pupil:', num2str(dDur), 'sec\n\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% alighment --- PROBABLY THERE ARE SOME ISSUES
%% front - scoring
PreX = ceil((vON/1000)/sBin);
Start = PreX * sBin; % start time
Score(1:PreX) = [];

%% front - pupil
preT = [1:length(Distance)]/fps + vON/1000;
idx = find(preT < Start);
Distance(idx) = [];
dMid(idx) = [];

sDur = length(Score)*sBin/60; % in min
pDur = (length(Distance)/fps+vON/1000)/60; % in min
dDur = (sDur - pDur)*60; % in sec
fprintf(['Score:', num2str(sDur), 'min\n']);
fprintf(['Distance:', num2str(pDur), 'min\n']);
fprintf(['Score - Pupil:', num2str(dDur), 'sec\n\n']);

%% alighment of end point
%% comparison
pEnd = length(Distance)/fps;
sEnd = sBin*length(Score);
postT = [1:length(Distance)]/fps;

if pEnd < sEnd % ephys data was longer than pupil monitoring    
    End = floor(pEnd/sBin);
    Score = Score(1:End);

    idx = find(postT >= length(Score)*sBin);
    Distance(idx) = [];
    dMid(idx) = [];
else % pupil monitoring was longer than ephys
    idx = find(postT >= sEnd);
    Distance(idx) = [];
    dMid(idx) = [];
end
    
sDur = length(Score)*sBin/60; % in min
pDur = (length(Distance)/fps)/60; % in min
dDur = (sDur - pDur)*60; % in sec
fprintf(['Score:', num2str(sDur), 'min\n']);
fprintf(['Distance:', num2str(pDur), 'min\n']);
fprintf(['Score - Pupil:', num2str(dDur), 'sec\n\n']);

%% time vectors in sec
vT = [1:length(Distance)]/fps; 
sT = sBin * [1:length(Score)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simple display
Yrange = [-3, 3];
cR = 0.9;
CO = [cR*ones(1,2) 1;
    cR 1 cR;
    1 cR*ones(1,2)]; % NREM, REM, AW
figure; subplot(2,4,[1:4]);
%% score
for s = 1:length(Score)
    if s == 1
        fill([0 sT(s) sT(s) 0], [Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)], 'k', 'FaceColor', CO(Score(s),:), 'LineStyle', 'none');hold on;
    else
        fill([sT(s-1) sT(s) sT(s) sT(s-1)], [Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)], 'k', 'FaceColor', CO(Score(s),:), 'LineStyle', 'none');
    end  
end

plot(vT, Distance, 'k-'); % diameter
% plot(vT, dMid, 'c-'); % eye center
hold off;
xlabel('sec');ylabel('z score');
box off;
ylim(Yrange);
xlim([0, sT(end)]);
title('pupil diameter');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% transitions analysis
%% var conversions
MyTs = vT;
Fs = fps;
NormFs = Distance;
Twin = sBin;

%%
tmp = diff(Score);
Tidx = find(tmp ~= 0); % when state changes?
idx = find(Tidx*Twin > MyTs(end)-Win);
Tidx(idx) = [];
nTs = length(Tidx); % # of transitions
nWins = 2*Win*Fs+1;
T = -1*Win:1/Fs:Win;

% matrices to build
Ttype = zeros(nTs,2); % transition type (before, after)
Tsignals = zeros(nTs,nWins); % pupil signals
StateSignals.sws = [];
StateSignals.rem = [];
StateSignals.aw = [];


for t = 1:nTs
    MyScore = Score(Tidx(t):Tidx(t)+1); % before, after
    MyTtimeON = Tidx(t) * Twin - Win;
    Idx = find(MyTs >= MyTtimeON, 1);
    MySignals = NormFs(Idx:Idx+nWins-1);
    
    % update
    Ttype(t,:) = MyScore;
    Tsignals(t,:) = MySignals;
    
    % signals
    if t > 1
        MyON = Tidx(t-1) * Twin;
        MyOFF = Tidx(t) * Twin;
        MyState = Score(Tidx(t));
        Idx = find(MyTs > MyON & MyTs < MyOFF);
        MySig = NormFs(Idx);
        if MyState == 1 % NREM
            StateSignals.sws = [StateSignals.sws; MySig];
        elseif MyState == 2 % REM
            StateSignals.rem = [StateSignals.rem; MySig];
        else
            StateSignals.aw = [StateSignals.aw; MySig];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% state-transitions
%% AW --> NREM (3 --> 1)
idx = find(Ttype(:,1) == 3 & Ttype(:,2) == 1);
if ~isempty(idx)
    My = Tsignals(idx,:);
    subplot(245);
    fill([-Win zeros(1,2) -Win],[Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)],'k','FaceColor', CO(3,:), 'LineStyle', 'none'); hold on;
    fill([0 Win*ones(1,2) 0],[Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)],'k','FaceColor', CO(1,:), 'LineStyle', 'none');    
    plot(T, My', 'Color', 0.85*ones(3,1));
    plot(T, mean(My),'k');
    plot(zeros(2,1),Yrange,'k:');hold off;
    axis([T(1) T(end) Yrange]);
    box off;
    title(['AW -> NREM: n=', num2str(length(idx))]);
    xlabel('sec');ylabel('normalized pupil diameter');
end

%% NREM --> REM (1 --> 2)
idx = find(Ttype(:,1) == 1 & Ttype(:,2) == 2);
if ~isempty(idx)
    My = Tsignals(idx,:);
    subplot(246);
    fill([-Win zeros(1,2) -Win],[Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)],'k','FaceColor', CO(1,:), 'LineStyle', 'none'); hold on;
    fill([0 Win*ones(1,2) 0],[Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)],'k','FaceColor', CO(2,:), 'LineStyle', 'none');
    plot(T, My', 'Color', 0.85*ones(3,1));
    plot(T, mean(My),'k');
    plot(zeros(2,1),Yrange,'k:');hold off;
    axis([T(1) T(end) Yrange]);
    box off;
    title(['NREM -> REM: n=', num2str(length(idx))]);
    xlabel('sec');
end
    
%% NREM --> AW (1 --> 3)
idx = find(Ttype(:,1) == 1 & Ttype(:,2) == 3);
if ~isempty(idx)
    My = Tsignals(idx,:);
    subplot(247);
    fill([-Win zeros(1,2) -Win],[Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)],'k','FaceColor', CO(1,:), 'LineStyle', 'none'); hold on;
    fill([0 Win*ones(1,2) 0],[Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)],'k','FaceColor', CO(3,:), 'LineStyle', 'none');
    plot(T, My', 'Color', 0.85*ones(3,1));
    plot(T, mean(My),'k');
    plot(zeros(2,1),Yrange,'k:');hold off;
    axis([T(1) T(end) Yrange]);
    box off;
    title(['NREM -> AW: n=', num2str(length(idx))]);
    xlabel('sec');
end

%% REM --> AW (2 --> 3)
idx = find(Ttype(:,1) == 2 & Ttype(:,2) == 3);
if ~isempty(idx)
    My = Tsignals(idx,:);
    subplot(248);
    fill([-Win zeros(1,2) -Win],[Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)],'k','FaceColor', CO(2,:), 'LineStyle', 'none'); hold on;
    fill([0 Win*ones(1,2) 0],[Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)],'k','FaceColor', CO(3,:), 'LineStyle', 'none');
    plot(T, My', 'Color', 0.85*ones(3,1));
    plot(T, mean(My),'k');
    plot(zeros(2,1),Yrange,'k:');hold off;
    axis([T(1) T(end) Yrange]);
    box off;
    title(['REM -> AW: n=', num2str(length(idx))]);
    xlabel('sec');
end


Ipath = 'I:\Science\SIPBS\sakata\Shuzo\Misc\Tomomi\Data\Ephys\P wave rec\with silicon probe\sponta pons rec\figs';
MySize = [18 8];
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', MySize); % width, height (in inches)
set(gcf, 'PaperPosition', [0 0.1 MySize-0.2]); %[left, bottom, width, height]
print('-djpeg',fullfile(Ipath,['Pupil_', Filebase,'_stateDependency.jpg']),'-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% overall signal stats
Sigs = [];
Grp = [];
MySig = StateSignals.sws;
Sigs = [Sigs; MySig];
Grp = [Grp; ones(length(MySig),1)];

MySig = StateSignals.rem;
Sigs = [Sigs; MySig];
Grp = [Grp; 2*ones(length(MySig),1)];

MySig = StateSignals.aw;
Sigs = [Sigs; MySig];
Grp = [Grp; 3*ones(length(MySig),1)];

figure;
subplot(121);
boxplot(Sigs,Grp,'notch','on');
box off;
set(gca,'XTickLabel',{'NREM','REM','AW'});
xlabel('states');ylabel('normalized pupil diameter');

CO = [0, 0, 1; % NREM (b)
    0, 1, 0; % REM (g)
    1, 0, 0]; % AW (r)

subplot(122);
x = Yrange(1):0.1:Yrange(2);
for s = 1:3
    N = hist(Sigs(Grp == s), x);
    % bar(x, N/sum(N), 'FaceColor', CO(s,:)); hold on;
    plot(x, N/sum(N), 'Color', CO(s,:)); hold on;
end
hold off;
box off;
xlabel('normalized pupil diameter');
ylabel('fraction');
legend({'NREM','REM','AW'});
xlim(Yrange);

MySize = [10 5];
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', MySize); % width, height (in inches)
set(gcf, 'PaperPosition', [0 0.1 MySize-0.2]); %[left, bottom, width, height]
print('-djpeg',fullfile(Ipath,['Pupil_', Filebase,'_stats.jpg']),'-r300');