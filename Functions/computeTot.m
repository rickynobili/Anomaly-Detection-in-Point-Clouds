function [tot] = computeTot(h_star,count,LCS,tot,index,sigma_pcd)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for h_s = 1:length(h_star)
    for a=2:count
        switch h_s
            case 1
                if (LCS(a,1)<=h_star(h_s) && LCS(a,1)>=0 && LCS(a,2) <= h_star(h_s) && LCS(a,2)>=0 && abs(LCS(a,3)) <= max(3*sigma_pcd,h_star(h_s))) 
                    if a<= index
                        tot(a-1) = tot(a-1)+1;
                    else
                        tot(a) = tot(a)+1;
                    end
                end
                
            case 2
                if (LCS(a,1)>=-h_star(h_s) && LCS(a,1)<=0 && LCS(a,2) <= h_star(h_s) && LCS(a,2)>=0 && abs(LCS(a,3)) <= max(3*sigma_pcd,h_star(h_s)))
                    if a<= index
                        tot(a-1) = tot(a-1)+1;
                    else
                        tot(a) = tot(a)+1;
                    end
                end
                
            case 3
                if (LCS(a,1)>=-h_star(h_s) && LCS(a,1)<=0 && LCS(a,2) >= -h_star(h_s) && LCS(a,2)<=0 && abs(LCS(a,3)) <= max(3*sigma_pcd,h_star(h_s)))
                    if a<= index
                        tot(a-1) = tot(a-1)+1;
                    else
                        tot(a) = tot(a)+1;
                    end
                end
        
            case 4
                if (LCS(a,1)<=h_star(h_s) && LCS(a,1)>=0 && LCS(a,2) >=-h_star(h_s) && LCS(a,2)<=0 && abs(LCS(a,3)) <= max(3*sigma_pcd,h_star(h_s))) 
                     if a<= index
                        tot(a-1) = tot(a-1)+1;
                    else
                        tot(a) = tot(a)+1;
                     end
                end
        end
    end
end

end

