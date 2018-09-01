function ksalt = kSalt(T)
        % T    Temperature, Celcius
        % m   Molality of the dissolved salt (mol/kg H2O)
        m = 4;
        T = T - 273.15;
        coef = [+0.11572, -6.0293e-4, +3.5817e-6, -3.7772e-9]; % Bakker (2003)      
%         T = T;
%         coef = [+0.257555, -0.157492e-3, -0.253024e-5, +0.438362e-8];
        temp = 0;
        for i = 1 : numel(coef)
            temp = temp + coef(i).*T.^(i-1);
        end
        ki = temp;
        ksalt = exp(ki.*m);
end

