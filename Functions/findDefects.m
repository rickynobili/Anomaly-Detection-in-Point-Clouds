function [found_defects] = findDefects(num,cmatrix,tmatrix,pc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global K;

found_defects = ones(1,length(pc.Location));

for a =1:length(cmatrix)
    
     if ((cmatrix(a,1) > (num))) % (floor(mean(tot)-std(tot)))/2)
         found_defects(a) = 0;
         if((cmatrix(a,1) == num+1)||(cmatrix(a,1) == num-1))
            cmatrix(a,:) = [0 0 1];
         else
             cmatrix(a,:) = [1 1 1];
         end
     end
end

% Se almeno il 75% dei punti vicini ad un punto considerato difetto non
% sono difetti allora quel punto non è un difetto
near_defects=0;

for x=1:pc.Count
    if(found_defects(x)==1)
        p = select(pc,x);
        KNNp = findNearestNeighbors(pc,p.Location,K);
        for y=1:K
            if(found_defects(KNNp(y))==1)
                near_defects = near_defects+1;
            end
        end
        if(near_defects<0.75*K)
            found_defects(x)=0;
        end
    end
end

ptc = pointCloud(pc.Location,'Color',cmatrix);
ptt = pointCloud(pc.Location,'Color',tmatrix);

subplot(1,2,1)
pcshow(ptt);
subplot(1,2,2)
pcshow(ptc);

end

