
% WRITE_VIDEO  Summary of this function goes here
% Author: Timothy Sipkens, 2020-04-03
%=========================================================================%

function [] = write_video(i,fname,cm,i_lim)

%-- Parse inputs --------------------%
if ~exist('cm','var'); cm = []; end
if isempty(cm); cm = gray; end

if ~exist('i_lim','var'); i_lim = []; end
%------------------------------------%

v = VideoWriter(fname);
v.FrameRate = 10;
open(v);

disp('Writing video file...');
tools.textbar(0);

for ii=1:1:size(i,3)
    imagesc(i(:,:,ii));
    axis image;
    colormap(cm);
    if ~isempty(i_lim); caxis(i_lim); end
    
    title(num2str(ii));
    drawnow;
    
    frames = getframe(gcf);
    
    if ~isempty(frames.cdata)
        writeVideo(v,frames);
    end
    tools.textbar(ii/size(i,3));
end
close(v);

end

