classdef SimpleWellComps < SimpleWell 

    methods
        function well = SimpleWellComps(W, varargin)
            % Class constructor
            well = well@SimpleWell(W, varargin{:});
        end
        %----------------------------------------------%
        function [weqs, ctrlEq, extra, extraNames, qMass, qVol, qm, wellSol] = computeWellEquationsComps(well, wellSol0, wellSol, resmodel, q_s, bh, packed, dt, iteration)
            % Compute well equations and well phase/pseudocomponent source terms
            [weqs, qMass, cqs, mix_s, status, cstatus, qVol, qm] = computeWellContributionsSingleWellComps(well, wellSol, resmodel, q_s, bh, packed);
            ctrlEq =  setupWellControlEquationsSingleWellComps(well, wellSol0, wellSol, bh, q_s, status, mix_s, resmodel);

            % Update well properties which are not primary variables
%             toDouble = @(x)cellfun(@double, x, 'UniformOutput', false);
%             cq_sDb = cell2mat(toDouble(qMass));

%             wellSol.cqs     = bsxfun(@rdivide, cq_sDb, resmodel.getSurfaceDensities);
            wellSol.cqs     = cqs;
            wellSol.cstatus = cstatus;
            wellSol.status  = status;
            extra = {};
            extraNames = {};
        end

    end
    
    
end