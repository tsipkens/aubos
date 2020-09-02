%=== generate_dots.m =====================================================%
% Function: Construct a dot pattern for the BOS background
% Author:   Samuel Grauer, 2017-10-25

function [varargout] = generate_dots(n,path,title,offset,params)
%-------------------------------------------------------------------------%
% Input:
%   w       Image width                                	[px]
%   h       Image height                               	[px]
%   n       Number of dots                            	[]
%   k       Color: w = white dots, k = black dots     	[]
%   path    Output path                                 []
%   title   Text                                        []
%   offset  Textbox offset                              [px]
%   params  Textbox parameters                          []
% Output:
%   I   Reference image                                 []
%-------------------------------------------------------------------------%


%--- Parse input ---------------------------------------------------------%
if ~exist('path','var'), f_save = 0;
else
    if isempty(path), f_save = 0;
    else, f_save = 1; end
end
if nargin < 3, f_title = 0; else, f_title = 1; end          % Title flag
if ~exist('offset','var'), offset = [25 25]; end            % Def. offset
if ~exist('params','var')                                   % Def. params
    params  = struct('boxcolor','w','boxopacity',0.825,...
        'fontsize',24,'font','Arial Bold');
end
%-------------------------------------------------------------------------%


%--- Construct figure ----------------------------------------------------%
f = figure;

set(f,'position',[50 50 1000 1000/sqrt(2)]);
set(gca,'position',[0 0 1 1]);
set(gca,'plotboxaspectratio',[1 1/sqrt(2) 1]);
set(gca,'visible','off');

pts = poissonDisc([1000 1000],5,n)';

fill(1000*[0 1 1 0],1000*[1 1 0 0],'w','linestyle','none'); hold on;
plot(pts(1,:),pts(2,:),'o','markerfacecolor','k');
axis off;
I = frame2im(getframe(gca));
close(f);

if f_title; I = insertText(I,offset,title,params); end
if f_save
    f = figure;
    imagesc(I);
    
    set(f,'position',[50 50 1189 1189/sqrt(2)]);
    set(gca,'position',[0 0 1 1]);
    set(gca,'plotboxaspectratio',[1 1 /sqrt(2) 1]);
    set(gca,'visible','off');
    
    set(f,'PaperUnits','centimeters');
    set(f,'PaperOrientation','portrait');
    set(f,'PaperPosition',[0 0 118.9 84.1]);
    ext = path(end-2:end);
    if strcmp(ext,'eps')
        print(f,path,'-deps','-r720');
    elseif strcmp(ext,'png')
        print(f,path,['-d',ext],'-r0');
    else
        error('Use png or pdf extension.');
    end
    close(f);
end
if nargout == 1; varargout{1} = im2double(I); end
%-------------------------------------------------------------------------%
end
%=== End =================================================================%






function [pts] = poissonDisc(sizeI,spacing,nPts,showIter)

% Purpose:
% N-dimensional poisson disc sampling function. This can also be used to
% randomly sample k pts from N-dimensional space with a minimum separation
% distance.
%
% Inputs:
% sizeI -   [required] Size of volume from which points are to be 
%           sampled
% spacing - [required] Minimum sepration distance between points
% nPts -    [Default is 0] if nPts = 0 For poisson disc sampling.
%           nPts = k to sample k-pts from N-dimensional space with
%           minimum separation distance
% showIter - [Default is 0] If showIter == 1, this option can be used to 
%            see how points are generated through each iteration. It can be
%            useful when code is taking a long time to generate points. 
%
% Output:
% pts - All eligible points
%
%
% Example:
% 1. Poisson disc sampling in 2-dimensional space. 
% sizeI = [512,512];
% spacing = 30;
% pts = poissonDisc(sizeI,spacing);
%
% 2. Sample k-pts in 3-dimensional space.
% sizeI = [512,512,192];
% spacing = 6;
% nPts = 10000;
% pts = poissonDisc(sizeI,spacing,nPts);
% 
% 3. Show iteration progress from poisson disc sampling in 2-dimension
% sizeI = [512,512];
% spacing = 6;
% nPts = 0;
% showIter = 1;
% pts = poissonDisc(sizeI,spacing,nPts,showIter);

% Mohak Patel, Brown University, 2016

%%%%%%% Initial parameters setup
% Parsing inputs and setting default values
if nargin == 3; showIter = 0; end
if nargin == 2; showIter = 0; nPts = 0; end

% Setting properties for iterations
ndim = length(sizeI);   % Number of Dimensions
k = 5;  % Number of 'dart' tries in each grid.
dartFactor = 4; %Select number of sample data in each iterations. Change it to
% reduce run time for code. Have to play around with number. 


%%%%%%% Making Grid read for iterations
%Make grid size such that there is just one pt in each grid
dm = spacing/sqrt(ndim);    % grize cell size [Bridson 2007]

%Make Grid
for i = 1:ndim
    sGrid{1,i} = 1:dm:sizeI(i);
end
[sGrid{:}] = ndgrid(sGrid{:});
sizeGrid = size(sGrid{1});

% Convert Grid points into a nx3 array;
for i = 1:ndim
    sGrid{i} = sGrid{i}(:);
end
sGrid = cell2mat(sGrid);

% Arrays to show eligible grids for dart throws and keeping score of darts
% thrown in a particular grid
emptyGrid = logical(ones(size(sGrid,1),1)); %Eligible Grids
nEmptyGrid = sum(emptyGrid);    %Number of eligible Grids
scoreGrid = zeros(size(emptyGrid)); %Score of darts thrown in Grid

% Darts to be thrown per iterations
% This hugely influences speed of the algorithm. Change dartFactor for it. 
if nPts == 0
    nPts = nEmptyGrid;
    ndarts = round(nEmptyGrid/dartFactor);
end
ndarts = round(nPts/dartFactor);

%%%%%%%%% Iterative process to generate points
% Initialize parameters
ptsCreated = 0;
pts = [];
iter = 0;

% Start Iterative process
tic
while ptsCreated<nPts & nEmptyGrid >0
    
    %Thrown darts in eligible grids
    availGrid = find(emptyGrid == 1);   %Eligible grids for dart throw
    dataPts = min([nEmptyGrid,ndarts]); % Darts to be thrown
    p = datasample(availGrid,dataPts,'Replace',false); %Select grids for darts
    tempPts = sGrid(p,:) + dm*rand(length(p),ndim); %Dart throw!!!
    
    
    % Find good dart throws
    [~,D] = knnsearch([pts;tempPts],tempPts,'k',2); %Finding distance between all darts(pts)
    D = D(:,2); 

    withinI = logical(prod(bsxfun(@lt,tempPts,sizeI),2)); %Eligible pts should be withing sizeI 
    eligiblePts = withinI & D>spacing; %elgible pts should also have minimum separation distance
    
    scorePts = tempPts(~eligiblePts,:); %Keep score from bad dart throws :(
    tempPts = tempPts(eligiblePts,:);   % Save good dart throws :)
    
    
    %Update empty Grid
    emptyPts = floor((tempPts+dm-1)/dm);
    emptyPts = num2cell(emptyPts,1);
    emptyIdx = sub2ind(sizeGrid,emptyPts{:});
    emptyGrid(emptyIdx) = 0;
    
    %Update score pts
    scorePts = floor((scorePts+dm-1)/dm);
    scorePts = num2cell(scorePts,1);
    scoreIdx = sub2ind(sizeGrid,scorePts{:});
    scoreGrid(scoreIdx) = scoreGrid(scoreIdx) + 1;
    
    %Update emptyGrid if scoreGrid has exceeded k dart throws
    emptyGrid = emptyGrid & (scoreGrid<k);
    
    %Update quantities for next iterations
    nEmptyGrid = sum(emptyGrid);
    pts = [pts;tempPts];
    ptsCreated = size(pts,1);
    iter = iter+1;
    ttoc = toc;
    
    %Display iteration details
    if showIter == 1
        disp(sprintf('Iteration: %d    Points Created: %d    EmptyGrid:%d    Total Time: %0.3f',iter,ptsCreated,nEmptyGrid,ttoc));
    end
    
end

% Cut down pts if more points are generated
if size(pts,1)>nPts
    p = 1:size(pts,1);
    p = datasample(p,nPts,'Replace',false);
    pts = pts(p,:);
end

end


