function [precision,recall,accuracy,F1,FPR] = ComputeMetrics(real_defect,found_defect)

TP = 0;
FN = 0;
FP = 0;
TN = 0;
count = length(real_defect);

for n = 1:count
    if (real_defect(n) == 1 && found_defect(n)==1)
        TP = TP + 1;
    end
    
    if (real_defect(n) == 1 && found_defect(n)==0)
        FN = FN + 1;
    end
    
    if (real_defect(n) == 0 && found_defect(n)==1)
        FP = FP + 1;
    end
    
    if (real_defect(n) == 0 && found_defect(n)==0)
        TN = TN + 1;
    end
end

precision = TP/(TP+FP);
recall = TP/(TP+FN);
accuracy = (TP+TN)/count;
F1 = 2*((precision*recall)/(precision+recall));
FPR = FP/(FP+TN);
end

