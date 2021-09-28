function [PHI] = computePHI(order,indici,LCS)

    PHI = [];
    phi = [];
    
    for i1=0:order
        i2=0;
        while i1+i2 <= order

            k=1:length(indici);

             phi = (LCS(indici(k),1).^i1).* (LCS(indici(k),2).^i2);

             

            if size(PHI)==0
                PHI = phi;
            else
                PHI = [PHI,phi];
            end

            i2 = i2+1;

        end
    end


end

