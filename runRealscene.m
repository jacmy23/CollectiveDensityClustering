clear;
clc;

curVideo = '.\example1';
curTrkName = 'klt_3000_10_trk_filt.txt';

curClipFrameSet = dir([curVideo '\*.jpg']);
curTrks = readTraks([curVideo '\' curTrkName]);

[XVset] = trk2XV(curTrks, 1, 2); % transform trajectories into point set

for ii = 2 : length(curClipFrameSet) - 2
    curFrame = imread([curVideo '\' curClipFrameSet(ii).name]);
    curFrame = im2double(curFrame);
    
    curIndex = find(XVset(:, 5) == ii);
    curX = XVset(curIndex,1:2);
    curV = XVset(curIndex,3:4);
    
    tic;
    
    x = curX;
    v = curV;
    
    index = find((v(:, 1).^2 + v(:, 2).^2).^ 0.5 == 0);
    x(index, :) = [];
    v(index, :) = [];
    
    N = size(x, 1);
    
    leg=sqrt(sum(v.^2, 2));
    vv=bsxfun(@rdivide,v,leg);
    dire = acos(vv(:, 1));
    dire(vv(:, 2) < 0) = bsxfun(@minus, 2 * pi,  dire(vv(:, 2) < 0));
    
    S = pdist2(x, x);
    Sd = pdist2(dire, dire);
    Sd(Sd > pi) = bsxfun(@minus, 2 * pi,  Sd(Sd > pi));
        
    xx = S(logical(triu(ones(N))));
    xx = xx(xx > 0);
    
    [NCLUST, cl, union_id] = CDC(N, S, Sd, xx, v);
    toc;
    
    if (NCLUST ~= 0)
        hold off;
        imshow(curFrame),hold on;
        
%         title(['Time = ' num2str(time_sum)], 'FontSize',20);
        
        ind = (cl == -1);
        scatter(x(ind,1),x(ind,2), 50, '+r');hold on;
        
        cluster_color = colormap(jet);
        colornum = length(unique(union_id));
        for k = 1 : NCLUST
            ic = int8((union_id(k) * 64) / (colornum * 1. + 2));
            id = find(cl == k);
            scatter(x(id, 1), x(id, 2), 50, cluster_color(ic, :), 'filled'); hold on;
            quiver(x(id, 1), x(id, 2),v(id, 1), v(id, 2), 1, 'color', cluster_color(ic, :), 'linewidth', 2.5);
        end
        drawnow;
%         pause;
    end
end




