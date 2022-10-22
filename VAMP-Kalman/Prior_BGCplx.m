classdef Prior_BGCplx
    
    properties
        rho
        gma
    end
    
    methods
        function obj = Prior_BGCplx(rho, gma)
            obj.rho = rho;
            obj.gma = gma;      
        end
          
        %{
        function s = generateRand(obj, N)
            sB = binornd(1, 1-obj.rho, [N 1]);
            sG = sqrt(1/(2*obj.gma)) * (randn([N 1]) + 1i*randn([N 1]));
            s = sB.* sG;
        end                
        %}
        
        function s = generateRand(obj, N)
            sB = zeros(N, 1);
            idx = randperm(N, ceil(N*(1-obj.rho)));
            sB(idx) = 1;
            sG = sqrt(1/(2*obj.gma)) * (randn([N 1]) + 1i*randn([N 1]));
            s = sB.* sG;
        end
        
        
        function [g, gPrime] = gfuncAndPrime(obj, r1, gma1)
            rhoBG = obj.rho;
            gmaBG = obj.gma;    
            
            gmaSP = gmaBG + gma1;            
            psi0 = (rhoBG/(1-rhoBG)) * (1+gma1/gmaBG) * exp(-gma1^2/gmaSP*abs(r1).^2);
                  
            moment1 = (gma1/gmaSP) * r1./ (1 + psi0);
            moment2 = (1/gmaSP) * (1 + (gma1^2/gmaSP)*abs(r1).^2)./ (1 + psi0);
            
            g = moment1;
            
            gPrime = gma1 * (moment2 - abs(g).^2);
        end  
    end
end

