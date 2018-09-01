nc = model.Mixture.numC;
time = zeros(numel(states),1);
nw = numel(wellSols{1});
qs = zeros(numel(states),nw);
sumqs = zeros(numel(states),1);
tempY = 0;
ztt = [];
acc =[];
for i = 1:numel(states)
    figure(1)
    p = plotCellData(model.G,convertTo(states{i}.pressure, barsa));
     plotWell(model.G, schedule.control(1).W);    
    title('Pressure, bar');
    xlabel('Pressure, bar');
    axis tight, view(35, 40), colorbar('SouthOutside');
    
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
        qOs(i, j) = convertTo(-ws(j).qOs , meter^3/day);
        qGs(i, j) = convertTo(-ws(j).qGs , meter^3/day);
         plot(time(1:i), qOs(1:i, j), lineSpec{:}, markerSpec{:});
          plot(time(1:i), qGs(1:i, j), lineSpec{:}, markerSpec{:});
          sw(:, i) = state.s(:, 1);
          so(:, i) = state.s(:, 2);
          st(:, i) = state.status(:);
          GOR(i, j) = ws(j).gor;
          if i == 14 || i == 19 || i ==23
              ztt = [ztt, state.zc(:, 1), state.zc(:, nc)];
              if isfield(state, 'ac')
                  acc = [acc, state.ac];
              end
          end
    end 
    qtt = [qOs(:, 1), qGs(:, 1)];
    % Axis
    xmax = convertTo(report.ReservoirTime(end), day);
    ymax = max(max(qGs(i,:)), tempY); 
    tempY = ymax;  
    axis([0 xmax 0 (ymax+10+0.1*ymax)]);
    title('Surface oil rate');
    xlabel('Time, days');
    ylabel('gas rate, m3/d');
    hold on
   pause(0.1);
end
