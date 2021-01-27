function TT_PupilVideoScreening(Filebase)

% 10 Oct 2018
% before DeepLabCut, this function will display which .avi file contains
% which segment. Thereby, you can choose an optimal an .avi file for
% training.

% NOTE: this function has not implemented a way to compensate frame dropping 
% but it is ok for a rough evaluation to determine an optimal video file

close all;

%% sleep scoring
load('manual sleep scoring\sleepscoring_results.mat');
sBin = 4; % in sec

%% video on
tmp = load([Filebase, '.video.evt']);
vON = tmp(1,1); % in ms

%% avi files
fps = 25;
vFolder = 'pupil video';
Vlist = ls(fullfile(vFolder, [Filebase, '*.avi']));
nFiles = size(Vlist,1);

Frames = zeros(nFiles,1);
for f = 1:nFiles
    MyFile = Vlist(f,:);
    fprintf([MyFile,'\n']);
    Vstr{f} = VideoReader(fullfile(vFolder,MyFile));
    Frames(f) = Vstr{f}.Duration*Vstr{f}.FrameRate;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% alighment
%% front - scoring
PreX = ceil(vON/1000/sBin);
Start = PreX * sBin; % start time
Score(1:PreX) = [];

%% front - pupil
preT = 1/fps*[1:Frames(1)] + vON/1000;
idx = find(preT < Start);
Frames(1) = Frames(1)-length(idx);

%% end
postT = 1/fps*[1:sum(Frames)]; 
End = floor(postT(end)/sBin);
Score = Score(1:End);

idx = find(postT >= End*sBin);
Frames(end) = Frames(end) - length(idx);

%% time vectors in sec
vT = 1/fps*[1:sum(Frames)]; 
sT = sBin * [1:length(Score)];

Grid = cumsum(Frames);
gT = 1/fps*Grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% display
%% simple display
Yrange = [0, 1];
cR = 0.9;
CO = [cR*ones(1,2) 1;
    cR 1 cR;
    1 cR*ones(1,2)]; % NREM, REM, AW
figure; 
%% score
for s = 1:length(Score)
    if s == 1
        fill([0 sT(s) sT(s) 0], [Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)], 'k', 'FaceColor', CO(Score(s),:), 'LineStyle', 'none');hold on;
    else
        fill([sT(s-1) sT(s) sT(s) sT(s-1)], [Yrange(1)*ones(1,2), Yrange(2)*ones(1,2)], 'k', 'FaceColor', CO(Score(s),:), 'LineStyle', 'none');
    end  
end

set(gca,'XTick', gT, 'XTickLabel', 1:length(gT));
hold off;grid on;
xlabel('file#');ylabel('z score');
box off;
ylim(Yrange);
xlim([0, sT(end)]);


%% image out
% MySize = [15 3];
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', MySize); % width, height (in inches)
% set(gcf, 'PaperPosition', [0 0.1 MySize-0.2]); %[left, bottom, width, height]
% print('-djpeg',fullfile('figs',[Filebase,'_VideoFileMapping.jpg']),'-r300');