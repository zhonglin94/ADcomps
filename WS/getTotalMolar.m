function state = getTotalMolar(model,state)

     numC = model.Mixture.numC;
     o = model.operators;
     rhom = state.rhom;
     s= state.s;
     ac = state.ac;
     xc= state.xc;
     yc = state.yc;
     
     
     Mw = o.pv.*rhom(:, 1).*s(:, 1).*(1-ac);
     Maco2 = o.pv.*rhom(:, 1).*s(:, 1).*ac;
     
     Mwt = sum(Mw);
     Maco2t = sum(Maco2);
     
     for i = 1 :  numC
       Mh(:, i) = o.pv.*rhom(:, 2).*s(:, 2).*xc(:, i) + o.pv.*rhom(:, 3).*s(:, 3).*yc(:, i);
       Mht(1, i) = sum(Mh(:, i));
     end
     Mt =  [Mwt, Mht];
     
     if isfield(model.inputdata, 'initMolar')
         iMt = model.inputdata.initMolar;
     else
         iMt = Mt;
     end
          
     RF = (iMt - Mt)./(iMt);
     
     state.Mt = Mt;
     state.RF = RF; 
     state.Mco2 = Maco2t;
end