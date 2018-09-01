%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
% Email: Yzhonglin_edu@163.com
% ADcomps is a simulator which is able to model the evaluate the development performance of CWI 
% through integrating appropriate equation of state and mixing into conventional compositional simulator. 
% Finding a creative algorithm that correlates phase behavior of oil-water system and solving non-linear 
% flow equations through Newton iteration efficiently are the two most challenging parts. The simulator
 % followed the coding convention of MRST which is a valuable open source simulator embedding in MATLAB. 
% This simulator is developed base on the frame of MRST 2017a and it should be run with the corresponding version.
%
%}

clear;
clc;
close all;
% We begin by loading the required modules
mrstModule clear
mrstModule add ADcomps-master ad-core ad-props mrst-gui  ad-blackoil
[ state, model, schedule] = compsExampleSetUp;
% Solve
nls = NonLinearSolver('useLinesearch', true);
[wellSols, states, report] = simulateScheduleAD(state, model, schedule, 'nonlinearsolver', nls);
% Plot after simulation
time = zeros(numel(states),1);
for i = 1:numel(states)
   
     figure(2)
     pp = plotCellData(model.G,states{i}.ac.*100);  
      plotWell(model.G, schedule.control(1).W);    
     title('Molar fraction, %');
    axis tight, view(35, 40), colorbar('SouthOutside');
    figure(3) % Plot for gas rate and cumulative gas versus time    
       % Line specification
    lineSpec = {'LineWidth', 2}; 
    markerSpec = {'Marker','o',  'MarkerFaceColor','green', 'MarkerEdgeColor', 'blue'};
    time(i) = convertTo(report.ReservoirTime(i), day);
    state =states{i};
    ws = wellSols{i};
    for  j = 1 : 1
        qOs(i, j) = convertTo(-ws(j).qOs , meter^3/day); % Surface oil rate
        qGs(i, j) = convertTo(-ws(j).qGs , meter^3/day); % Surface Gas rate
        qWs(i, j) = convertTo(-ws(j).qWs , meter^3/day); % Surface water rate
        
        plot(time(1:i), qOs(1:i, j), lineSpec{:}, markerSpec{:}); % Plot surface oil rate
    end 
    % Axis
    xmax = convertTo(report.ReservoirTime(end), day);
    ymax = max(qOs(i,:)); 
    tempY = ymax;  
    axis([0 xmax 0 600]);
    title('Surface oil rate');
    xlabel('Time, days');
    ylabel('oil rate, m3/d');
    hold on
   pause(0.05);
end

