function [OUT, params, H5list] = DLC_PupilReconstruction(Filebase)

% 10 Oct 2018
% after analyzing video files to get .h5 files, this function integrate all
% coordination information

% OUT.mat (time x elements)
%   .L time x coordinate(x,y)
%   .R 
%   .D distance (pixels)
%   .t (in sec)
%
% params.coorpos ... position extracted from .mat
% params.fps ... frames/sec
%

%% specific params 
CoorPosition = [1, 2, 4, 5]; % L(x,y), R(x,y)
fps = 25;
params.coorpos = CoorPosition;
params.fps = fps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% frame drop analysis
Vlist = ls([Filebase, '*.avi']);
nFiles = size(Vlist,1);

%% checking duration
Durs = zeros(nFiles,1);
for f = 1:nFiles
    MyFile = Vlist(f,:);
    vObj = VideoReader(MyFile);    
    Durs(f) = vObj.Duration;
end

MaxDur = max(Durs);
DroppedFrames = zeros(nFiles,1);
for f = 1:nFiles-1 % ignore the final file
    MyFile = Vlist(f,:);
    vObj = VideoReader(MyFile);
    
    dDur = MaxDur - vObj.Duration;
    DroppedFrames(f) = dDur*fps;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% h5 file list
H5list = ls([Filebase, '*.h5']);
nFiles = size(H5list,1);

%% construct coordination matrix
BigMat = [];
for f = 1:nFiles
    MyFile = H5list(f,:);
    fprintf([MyFile,'\n']);
    data = h5read(MyFile,'/df_with_missing/table');
    MyInfo = data.values_block_0(CoorPosition,:);
    
    %% insertion missed frame (randomly insert missed frames)
    if DroppedFrames(f) ~= 0
        nFrames = size(MyInfo,2);
        if nFrames > 1 % in rare case, there is a single frame in a video file, which causes an error
            Insert = floor(nFrames/(DroppedFrames(f)+1));
            Idx = 0:Insert:nFrames;
            Idx(1) = []; Idx(end) = [];
            for d = 1:DroppedFrames(f)
                Pre = MyInfo(:,Idx(d)+(d-1)-1);
                Post = MyInfo(:,Idx(d)+(d-1));
                My = mean([Pre, Post],2);
                MyInfo = [MyInfo(:,1:Idx(d)+(d-1)-1), My, MyInfo(:,Idx(d)+(d-1):end)];
            end
        end
    end
 
    %% update
    BigMat=[BigMat, MyInfo];
end

%% extract each coordinate and distance
L = BigMat(1:2,:)';
R = BigMat(3:4,:)';
D = sqrt((L(:,1)-R(:,1)).^2+(L(:,2)-R(:,2)).^2);

%% visualization
t=(1/fps)*[1:1:length(D)]; % in sec
plot(t/60,D)
xlabel('min');ylabel('pupil diameter (pixel)');
box off;

%% output
OUT.mat = BigMat';
OUT.L = L;
OUT.R = R;
OUT.D = D;
OUT.t = t;

save([Filebase,'.pupil.mat'], 'OUT', 'params', 'H5list','DroppedFrames');
