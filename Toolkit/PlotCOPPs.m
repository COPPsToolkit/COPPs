function [hFig] = PlotCOPPs(projections, params)

% Check if plotting is requested
if ~isfield(params.plotting,'NoPlot')
elseif isfield(params.plotting,'NoPlot') & params.plotting.NoPlot
    return;
end

%% Set options

if ~isfield(params,'plotting')
params.plotting.Default = 1;
end

if ~isfield(params.plotting,'first_date')
    nodates = 1;
else
    nodates = 0;
    if ~isfield(params.plotting,'freq')
        warning('Frequency of data must be given: month/quarter/year')
    end
    params.plotting.dates = dateshift(params.plotting.first_date,'start',params.plotting.freq,0:params.PastPeriods+params.T_full);
end

if ~isfield(params.plotting,'P_past')
    P_past   = 20;
else
    P_past = params.plotting.P_past;
end
if ~isfield(params.plotting,'P_future')
    P_future = 36;
else
    P_future = params.plotting.P_future;
end

if ~isfield(params.plotting,'Data')
    Data = {projections.simple_rules.data  , '-'  , 1 ,[0,0,0], 'Baseline' ; ...
        projections.COPPs.data             , '--' , 1 ,[0,0,1], 'COPPs'; ...
        };
else
    Data = params.plotting.Data;
end

if ~isfield(params.plotting,'VarsToPlot')
    for i_V2P = 1:length({params.PolicyInstrumentsAndShocks{:,1}})
        VarsToPlot{i_V2P,1} = params.PolicyInstrumentsAndShocks{i_V2P,1};
    end
    for i_V2P = 1:length(params.LossFunctionVariables)
        VarsToPlot{length({params.PolicyInstrumentsAndShocks{:,1}})+i_V2P,1} = params.LossFunctionVariables{i_V2P,1};
    end
else
    VarsToPlot = params.plotting.VarsToPlot;
end

if nodates == 1
    if ~isfield(params.plotting,'VerticalLines')
        VerticalLines = 0;
    else
        VerticalLines = params.plotting.VerticalLines;
    end
elseif nodates == 0
    if ~isfield(params.plotting,'VerticalLines')
        VerticalLines = params.plotting.dates(1,params.PastPeriods);
    else
        VerticalLines = params.plotting.dates(1,VerticalLines);  % these need to be in a column vector
    end
end

%% Create plots
if nodates~=1
    dates = params.plotting.dates(params.PastPeriods-P_past+1:params.PastPeriods + P_future); 
end
nPlotLines = size(Data,1);
nPlotVars  = size(VarsToPlot,1);
PlotVarIndexes = nan(1,nPlotVars);

for i = 1:nPlotVars
    PlotVarIndexes(i) = find(strcmp(projections.simple_rules.model.VarNames,VarsToPlot{i,1}));
end

if nPlotVars == 1
    mp = 1;
    np = 1;
elseif nPlotVars <= 4
    mp = 2;
    np = 2;
elseif nPlotVars <= 6
    mp = 2;
    np = 3;
else
    mp = 3;
    np = 3;
end

hAxis = cell(nPlotVars,1);

for i = 1:nPlotVars
    
    subplotidx = mod(i,9);
    if subplotidx == 0
        subplotidx = 9;
    end
    if subplotidx == 1
        width  = 12*np;
        height = 8*mp;
        hFig = figure('Color',            [1 1 1], ...
            'Units',            'centimeters', ...
            'PaperUnits',       'centimeters', ...
            'PaperSize',        [width height], ...
            'PaperPosition',    [-0.5 , 0.5 , width-1 , height-1], ...
            'NumberTitle',      'off', ...
            'Visible',          'on', ...
            'Name',             'COPPs');
        hold on;
    end
    
    subplot(mp,np,subplotidx);
    hAxis{i} = gca;
    hold on;
    
    varYLim = [nan nan];
    for j = 1:nPlotLines
        lineData          = Data{j,1};
        lineStyle         = Data{j,2};
        
        lineDeviations = [lineData.Past.Deviations(PlotVarIndexes(i),(end-P_past+1):end) , lineData.Forecast.Deviations(PlotVarIndexes(i),1:P_future)];
        lineTrends     = [lineData.Past.Trends(PlotVarIndexes(i),(end-P_past+1):end)     , lineData.Forecast.Trends(PlotVarIndexes(i),1:P_future)];
        lineLevels     = lineTrends + lineDeviations;
        
        if nodates==1
            hLines(j) = plot((1-P_past):P_future, lineLevels, lineStyle ,'LineWidth',Data{j,3},'Color',Data{j,4});
        else
            hLines(j) = plot(dates, lineLevels', lineStyle ,'LineWidth',Data{j,3},'Color',Data{j,4});

            if params.plotting.freq            == 'quarter'
                hAxis{i}.XAxis.TickLabelFormat ='yyyy';
            elseif params.plotting.freq        == 'month'
                hAxis{i}.XAxis.TickLabelFormat ='mm-yyyy';
            elseif params.plotting.freq        == 'year'
                hAxis{i}.XAxis.TickLabelFormat ='yyyy';
            end
        end


        if isfield(params.plotting,'VarsToPlot')
            varYLim  = params.plotting.VarsToPlot{i,2};
        else
            varYLim = [min([varYLim(1),lineLevels]) max([varYLim(2),lineLevels])];
        end
        
    end
    
    box on;
    
    
    if isfield(params.plotting,'VarsToPlot')
        varLabel = params.plotting.VarsToPlot{i,3};
    else
        varLabel = VarsToPlot(i);
    end
    title(varLabel,'Interpreter','none', 'FontSize', 12 , 'FontWeight','Normal');
    ylim(varYLim);
    
        
    for vl = 1:length(VerticalLines)
        xdata = [VerticalLines(vl) , VerticalLines(vl)];
        ydata = varYLim;
        plot(xdata,ydata,'-','Color',[0.5,0.5,0.5],'LineWidth',0.5);
    end
    
    if nodates==1
        plot((1-P_past):P_future, zeros(P_past+P_future) , 'Color',[0,0,0,.5])
        plot((1-P_past):P_future, lineTrends ,'LineStyle',':', 'Color',[0,0,0,.5])
    else
        plot(dates, zeros(1,length(dates)) , 'Color',[0,0,0,.5])
        plot(dates, lineTrends ,'LineStyle',':', 'Color',[0,0,0,.5])
    end
end
hLegend = legend(hLines,Data{:,5},'Location','Best');

end