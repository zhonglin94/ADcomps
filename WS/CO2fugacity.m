function fug = CO2fugacity(model, ac, P)
% Get CO2 fugacity in aqueous phase

% PR with WS mixing rules requires three binary interaction parameters:
%the binary interaction parameter k, and NRTL parameters
%tau (interaction energy parameter) and alpha (non-randomness parameter)
%All of them can be temperature-dependant
% del.k1  = [0 0.3073;0.3073 0]; %Binary interaction parameter k
% del.k2 = [0 0.1141;0.1141 0]; %NRTL parameter alpha
% del.k3 = [0 4.3870;0.3930 0]; %NRTL paramter tau
del.k1 = [0 0.3073;0.3073 0];
del.k2 = [0 0.1120;0.1120 0];
del.k3 = [0 4.3570;0.4130 0];

cCO2 =ac;
cH2O = 1-ac;
comp = {cCO2, cH2O};
id = [3, 61];
mix = getCompsProps(id);
mix.del = del;
T = model.T;
fug = fugacityPRWS(comp, mix, P, T); 
end

function fug = fugacityPRWS(comp, mix, P, T)
R = 8.31; %Ideal gas constant, J/mol K
x = comp;
numC = numel(comp);
ng = numel(double(P));
Pc =   mix.Pc;      % Critical Pressures ( Pa )
Tc =   mix.Tc;          % Critical Temperatures ( K )
w  =   mix.w;         % Acentric Factors
del =   mix.del; % Interaction Coefficients 
k1 = del.k1;
k2 = del.k2;
k3 = del.k3;
%%
for i = 1:numC
   alpha(i) = alpha_function_ws(w(i), Tc(i), T); %Evaluates the alpha function (may differ in modifications of PR-EoS)
   a(i) = 0.45724*(R*Tc(i))^2/Pc(i)*alpha(i);
   b(i) = 0.0778*R*Tc(i)/Pc(i);
end

%Activity coefficient model
A_ex = 0;
g = zeros(numC,numC);
for i = 1:numC
    for j = 1:numC
        g(i,j) = exp(-k2(i,j).*k3(i,j));        
    end
end
for i = 1:numC
    Num_A_ex = 0;
    Den_A_ex = 0;
    for j = 1:numC
        Num_A_ex = Num_A_ex + x{j}.*k3(j,i).*g(j,i);
    end
    for k = 1:numC
        Den_A_ex = Den_A_ex + x{k}.*g(k,i);    
    end
    A_ex = A_ex + R.*T.*x{i}.*Num_A_ex./Den_A_ex;
end

%Mixing rules
C = log(sqrt(2) - 1)./sqrt(2);
D = A_ex./(C.*R.*T);
Q = 0;
b_mix = zeros(numC, numC);
for i = 1:numC
   D = D + x{i}.*a(i)./(b(i).*R.*T);
   for j = 1:numC
        b_mix(i,j) = ((b(i) - a(i)./(R.*T)) + (b(j) - a(j)./(R.*T)))./2.*(1-k1(i,j));
        Q = Q + x{i}.*x{j}.*b_mix(i,j);
   end
end
%Mixture parameters
am = R.*T.*Q.*D./(1-D);
bm = Q./(1-D);
A = am.*P./(R.*T)^2;
B = bm.*P./(R.*T);
%%
 % Get the roots and find an appropriate one for the specified phase
        phase =1;
        r = getCubicRoots(A,B);
        rootForSort = cell(1,3);
        for rootIndex = 1:3
               isImage = abs(imag(double(r{rootIndex})))>eps;
               r{rootIndex}(isImage) = NaN;
               isNegative = double(r{rootIndex})<-eps;
               r{rootIndex}(isNegative) = NaN;         
               rootForSort{rootIndex} = double(r{rootIndex});        
        end

        [~,minIndex] = min(cell2mat(rootForSort),[],2);
        [~,maxIndex] = max(cell2mat(rootForSort),[],2);    
        % Phase = 0 --> GAS Phase,   phase = 1 --> LIQUID Phase
        Z = cell(ng,1);
        for i =1:ng
        Z{i}=phase.*r{minIndex(i)}(i)+(1-phase).*r{maxIndex(i)}(i); 
            if ~isreal(double(Z{i}))
                % Keeping Z a real number in case the image part of Z is
                % smaller than eps
                if isa(Z{i}, 'ADI')
                    % handle ADI class
                    Z{i}.val = real(Z{i}.val);
                    for j = 1:numel(Z{i}.jac)
                        Z{i}.jac(i) = num2cell(real(cell2mat(Z{i}.jac(j))), [1 2]);
                    end
                else
                    Z{i} = real(Z{i});
                end
            end
        end
        Z = vertcat(Z{:});
%%
%.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*
%Calculates the fugacity coefficient
%.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*.*
%Molar volume
v = Z.*R.*T./P;

% f = zeros(numC,1);
[fug, f] = deal(cell(1,numC));
for comp = 1:numC
    %Activity coefficient
    Num_A_ex = 0;
    Den_A_ex = 0;
    for j = 1:numC
        Num_A_ex = Num_A_ex + x{j}.*k3(j,comp).*g(j,comp);
    end
    for k = 1:numC
        Den_A_ex = Den_A_ex + x{k}.*g(k,comp);    
    end
    term_1 = Num_A_ex./Den_A_ex;

    term_2 = 0;
    for j = 1:numC
        Num_A_ex = 0;
        for l = 1:numC
           Num_A_ex = Num_A_ex + x{l}.*k3(l,j).*g(l,j);
        end
        Den_A_ex = 0;
        for k = 1:numC
            Den_A_ex = Den_A_ex + x{k}.*g(k,j);  
        end
        term_2 = term_2 + x{j}.*g(comp,j)./Den_A_ex.*(k3(comp,j) - Num_A_ex./Den_A_ex);
    end
    log_gamma = term_1 + term_2;
    
    %Auxiliary derivatives
    der2Q_dn2 = 0;
    for j = 1:numC
        der2Q_dn2 = der2Q_dn2 + 2.*x{j}.*b_mix(comp,j);
    end

    derD_dn = a(comp)./(b(comp).*R.*T) + log_gamma./C;

    derb_dn = 1./(1-D).*der2Q_dn2 - Q./(1-D).^2.*(1-derD_dn);

    der2a_dn2 = R.*T.*(D.*derb_dn + bm.*derD_dn);

    %Fugacity coefficient
    term1 = -log(P.*(v-bm)./(R.*T));
    term2 = 1./bm.*derb_dn.*(Z-1);
    term3 = 1./(2.*sqrt(2)).*am./(bm.*R.*T);
    term4 = (1./am.*der2a_dn2 - 1./bm.*derb_dn);
    term5 = log((v+bm.*(1-sqrt(2)))./(v+bm.*(1+sqrt(2))));
    f{comp} = exp(term1 + term2 + term3.*term4.*term5);
    fug{comp} = x{comp}.*P.*f{comp};
end
end 

 function roots = getCubicRoots(A,B)
        a =  1;
        b = -(1-B);
        c =  A-3.*B.^2-2.*B;
        d = -(A.*B-B.^2-B.^3);
        x1 = (((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6*a.^2)).^2 + (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)... 
        - b./(3.*a) - (- b.^2./(9.*a.^2) + c./(3.*a))./(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6*a.^2)).^2 + (c./(3.*a) - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3)...
         - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3);
        x2 =  (- b.^2./(9.*a.^2) + c./(3.*a))./(2.*(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (c./(3.*a) - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3)...
         - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)) + (3.^(1./2).*((- b.^2./(9.*a.^2) + c./(3.*a))./(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (c./(3.*a)...
         - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3) + (((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + ...
        (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)).*i)./2 - b./(3.*a) - ...
        (((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)./2;

        x3 = (- b.^2./(9.*a.^2) + c./(3.*a))./(2.*(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (c./(3.*a) - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3)... 
        - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)) - (3.^(1./2).*((- b.^2./(9.*a.^2) + c./(3.*a))./(((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 + (c./(3.*a)...
         - b.^2./(9.*a.^2)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3) + (((d./(2.*a) + b.^3./(27.*a.^3) - (b.*c)./(6.*a.^2)).^2 +... 
        (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)).*i)./2 - b./(3.*a) - (((d./(2.*a) + b.^3./(27.*a.^3) ...
        - (b.*c)./(6.*a.^2)).^2 + (- b.^2./(9.*a.^2) + c./(3.*a)).^3).^(1./2) - b.^3./(27.*a.^3) - d./(2.*a) + (b.*c)./(6.*a.^2)).^(1./3)./2;
        roots = {x1,x2,x3};
   end


 

