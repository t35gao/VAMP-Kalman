classdef Prior_GaussianCplx
    
    properties
        mu
        gma
    end
    
    methods
        function obj = Prior_GaussianCplx(mu, gma)
            obj.mu = mu;
            obj.gma = gma;      
        end
          
        
        function s = generateRand(obj, N)
            s = sqrt(1/2) * (randn([N 1]) + 1i*randn([N 1]));
            s = obj.mu + sqrt(1/obj.gma)*s;
        end                
        
        
        function [g, gPrime] = gfuncAndPrime(obj, r1, gma1)
            g = (obj.gma*obj.mu + gma1*r1) / (obj.gma + gma1);
            gPrime = gma1 / (obj.gma + gma1) * ones(size(r1));
        end  
    end
end

