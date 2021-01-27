function [Vstr, vMat] = DLC_MakeKwikVideo(Filebase, fps, Speed, nParts)

% 10 Oct 2018
% this function creats a down-sampled movie based on DeepLabCut results
%
% after making h5 files, run this function

% Speed ... how many times faster? (this function will down-sample accordingly)
% nParts ... how many body parts were labelled?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load h5 files
H5list = ls([Filebase, '*.h5']);
nFiles = size(H5list,1);

%% construct coordination matrix
fprintf('reading h5 files ... \n');
BigMat = [];
for f = 1:nFiles
    MyFile = H5list(f,:);
    fprintf([MyFile,'\n']);
    data = h5read(MyFile,'/df_with_missing/table');
    tmp = data.values_block_0;
    tmp(3:3:nParts*3,:) = []; % deleting likelihood values
    BigMat=[BigMat, tmp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% video list
fprintf('reading video files ... \n');
% Vlist = ls([Filebase, '*.wmv']); % .... WHEN YOU DEAL WITH WMV FILES
Vlist = ls([Filebase, '*.avi']);
nFiles = size(Vlist,1);

%% constructing video object structure
Vstr = [];
vMat = []; % elements x frames ... (1, video#; 2, frame#; 3, timeInFile)
Current = 0;
for f = 1:nFiles
    MyFile = Vlist(f,:);
    fprintf([MyFile,':']);
    
    %% to load video object
    Vstr{f} = VideoReader(MyFile);
    nFrames = floor(Vstr{f}.Duration*Vstr{f}.FrameRate);
    fprintf([num2str(nFrames), ' frames\n']);

    %% parameters (video file#, frame#, time)
    tmp = [f*ones(1,nFrames);
        [1:nFrames]+Current;
        [0:nFrames-1]/fps];
    
    %% update
    vMat = [vMat, tmp];
    Current = tmp(2,end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('creating a video file w/ markers ... \n');

%% sampling position
Idx = 1:Speed:vMat(2,end); % key for down-sampling
CO = round(hsv(nParts)*256); % colors for body parts
mSize = 1; % 1 ... add one around the center (if marks are too small, increase this
COmat = zeros(2*mSize+1,2*mSize+1,3,nParts); % where subset to replace colors
for p = 1:nParts
    for d = 1:3
        COmat(:,:,d,p) = repmat(CO(p,d), 2*mSize+1,2*mSize+1,1,1);
    end
end

%% make kwik movie
Vout = VideoWriter([Filebase, '_x', num2str(Speed),'.mp4'],'MPEG-4');
Vout.FrameRate = fps;
open(Vout);
for i = 1:length(Idx)
    if rem(i,fps*30) == 0
        fprintf('.');
    end
    MyIdx = Idx(i);
    MyFile = vMat(1, MyIdx); % file#
    % MyFrame = vMat(2, MyIdx);
    MyTime = vMat(3, MyIdx); % time
    
    %% extract frame
    Vstr{MyFile}.CurrentTime = MyTime;
    MyFrame = readFrame(Vstr{MyFile});
    
    %% marking parts
    MyCoo = floor(BigMat(:,MyIdx));
    for p = 1:nParts
        x = MyCoo(2*p);
        y = MyCoo(2*p-1);
        MyFrame(x-mSize:x+mSize,y-mSize:y+mSize,:) = squeeze(COmat(:,:,:,p));
    end
    
    %% update video
    writeVideo(Vout, MyFrame);
end
fprintf('\n');

close(Vout);

