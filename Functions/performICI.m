function [h_star] = performICI(h_star,sigma_teta_h,gammaICI,h1,quadrante,z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

D = zeros(length(h1),2);
        
    % Compute the intervals for the ICI
    for indice_z = 1:length(h1)
        D(indice_z,1) = z(indice_z)-gammaICI*sigma_teta_h;
        D(indice_z,2) = z(indice_z)+gammaICI*sigma_teta_h;
    end

    inferiore = D(1,1);
    superiore = D(1,2);

    for indice_d=1:length(h1)-1

        ni = D(indice_d +1 , 1);
        ns= D(indice_d +1 , 2);

        if (inferiore < ni && ni < superiore && superiore < ns)
            inferiore = ni;
            h_star(quadrante) = h1(indice_d+1);
        end

        if (superiore > ns && ni < inferiore && inferiore < ns)
            superiore = ns;
            h_star(quadrante) = h1(indice_d+1);
        end

        if (inferiore > ns) || (superiore < ni )
            h_star(quadrante) = h1(indice_d);
            break;
        end

        if (inferiore < ni) && (superiore > ns )
            inferiore = ni;
            superiore = ns;
            h_star(quadrante) = h1(indice_d+1);
        end

        if(indice_d == length(h1)-1)
             h_star(quadrante) = h1(indice_d+1);
             break;
        end
    end
end

