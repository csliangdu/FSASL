function [output] = prettyPlot(xData,yData,options)
% prettyPlot(xData,yData,options)
%
% available options:
% legend - cell array containing legend entries (default = [])
% title - string containing plot title (default = [])
% xlabel - string containing x-axis label (default = [])
% ylabel - string containing y-axis label (default = [])
% lineWidth - width of lines (default = 3)
% colors - cell array or (n by 3) matrix containing line colors (default = 'b')
% lineStyles - cell array containing line styles (default = '-')
% markerSize - size of markers (default = 12)
% markers - cell array containing markers (default = [])
% markerSpacing - (n by 2) matrix containing spacing between markers and offset for first marker
% xlimits - 2-vector containing lower and upper limits for x-axis (can be inf to show full range)
% ylimits - 2-vector containing lower and upper limits for y-axis (can be inf to show full range)
% logScale - can be 0 (regular scale), 1 (semilogx), 2 (semilogy), or 3 (loglog)
% legendLoc - location of legend (default = 'Best')
% useLines - whether to use lines (default = 1)
% fillFace - whether to fill markers (default = 1)
% errors - (n by p by 2) array containing upper and lower error lines
% errorStyle - line style for error bars
% errorColors - colors for error bars
% errorFill - set to 1 to use a filled-in error area instead of lines (default = 0)
% labelLiens - place legend values next to lines (default = 0)
% labelLoc - can be 0 (randomly placed), 1 (random, then furthest from
% existing), or 2 (furthest from other lines)
% labelRotate - rotates line labels (default = 1)

if nargin < 3
    options = [];
end

[legendStr,plotTitle,plotXlabel,plotYlabel,lineWidth,colors,lineStyles,...
    markerSize,markers,markerSpacing,xlimits,ylimits,logScale,legendLoc,...
    useLines,fillFace,errors,errorStyle,errorColors,errorFill,labelLines,labelLoc,labelRotate] = ...
    prettyPlotProcessOptions(options,'legend',[],'title',[],'xlabel',[],'ylabel',[],...
    'lineWidth',3,'colors',[],'lineStyles',[],...
    'markerSize',12,'markers',[],'markerSpacing',[],...
    'xlimits',[],'ylimits',[],...
    'logScale',0,'legendLoc','Best','useLines',1,'fillFace',1,...
    'errors',[],'errorStyle',{'--'},'errorColors',[],'errorFill',0,'labelLines',0,'labelLoc',2,'labelRotate',1);

if nargin < 2
    fprintf('\nUsage: prettyPlot(xData,yData,options)\n\nCalling with no arguments returns the default options\nUse the help function for details\n');
    options.legend = legendStr;
    options.title = plotTitle;
    options.xlabel = plotXlabel;
    options.ylabel = plotYlabel;
    options.lineWidth = 3;
    options.colors = colors;
    options.lineStyles = lineStyles;
    options.markerSize = markerSize;
    options.xlimits = xlimits;
    options.ylimits = ylimits;
    options.logScale = logScale;
    options.legendLoc = legendLoc;
    options.useLines = useLines;
    options.fillFace = fillFace;
    options.errors = errors;
    options.errorStyle = errorStyle;
    options.errorColors = errorColors;
    options.errorFill = 0;
    options.labelLines = 0;
    options.labelLoc = 2;
    options.labelRotate = 1;
    output=options;
    return
else
    output = [];
end

if logScale == 1
    plotFunc = @semilogx;
elseif logScale == 2
    plotFunc = @semilogy;
elseif logScale == 3
    plotFunc = @loglog;
else
    plotFunc = @plot;
end

if useLines == 0
    defaultStyle = 'b.';
else
    defaultStyle = 'b';
end

if iscell(yData)
    nLines = length(yData);
else
    nLines = size(yData,1);
end

for i = 1:nLines
    
    % Get yData for line
    if iscell(yData)
        y{i} = yData{i};
    else
        y{i} = yData(i,:);
    end
    
    % Get xData for line
    if isempty(xData)
        x{i} = 1:length(y);
    elseif iscell(xData)
        if isempty(xData{i})
            x{i} = 1:length(y{i});
        else
            x{i} = xData{i};
        end
    elseif size(xData,1) == 1
        x{i} = xData(1:length(y{i}));
    else
        x{i} = xData(i,:);
    end
    
    % Plot
    h(i) = plotFunc(x{i},y{i},defaultStyle);
    hold on;
end

if ~isempty(errors)
    if isempty(errorColors)
        errorColors = colors+.75;
        errorColors(:) = min(errorColors(:),1);
    end
    for i = 1:nLines
        if errorFill == 0
            hL = plotFunc(x{i},errors{i,1},defaultStyle);
            hU = plotFunc(x{i},errors{i,2},defaultStyle);
            applyStyle(hL,i,lineWidth,errorColors,errorStyle,[],markerSpacing);
            applyStyle(hU,i,lineWidth,errorColors,errorStyle,[],markerSpacing);
        else
            hE = fill([x{i}(:);flipud(x{i}(:))],[errors{i,1}(:);flipud(errors{i,2}(:))],errorColors(i,:),'LineStyle','none');
        end
    end
    
    % Re-plot the lines
    for i = 1:nLines
        h(i) = plotFunc(x{i},y{i},defaultStyle);
    end
end



if isempty(markerSpacing)
    for i = 1:nLines
        h(i) = applyStyle(h(i),i,lineWidth,colors,lineStyles,markers,markerSpacing);
    end
else
    for i = 1:nLines
        h(i) = applyStyle(h(i),i,lineWidth,colors,lineStyles,markers,markerSpacing);
        if ~isempty(markers) && ~isempty(markers{1+mod(i-1,length(markers))})
            hM = plotFunc(x{i}(markerSpacing(i,2):markerSpacing(i,1):end),y{i}(markerSpacing(i,2):markerSpacing(i,1):end),'b.');
            markerLocations{i} = markerSpacing(i,2):markerSpacing(i,1):length(x{i});
            applyStyle(hM,i,lineWidth,colors,[],markers,[]);
            hM = plotFunc(x{i}(markerSpacing(i,2)),y{i}(markerSpacing(i,2)),defaultStyle);
            h(i) = applyStyle(hM,i,lineWidth,colors,lineStyles,markers,[]);
        end
    end
end




set(gca,'FontName','AvantGarde','FontWeight','normal','FontSize',12);

if ~isempty(plotTitle)
    hTitle = title(plotTitle);
    set(hTitle,'FontName','AvantGarde','FontSize',12,'FontWeight','bold');
end

if ~isempty(plotXlabel) || ~isempty(plotYlabel)
    h1 = xlabel(plotXlabel);
    h2 = ylabel(plotYlabel);
    set([h1 h2],'FontName','AvantGarde','FontSize',12,'FontWeight','normal');
end

set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'LineWidth'   , 1         );

if ~isempty(xlimits)
    xl = xlim;
    xlimits(xlimits == -inf) = xl(1);
    xlimits(xlimits == inf) = xl(2);
    xlim(xlimits);
else
    xlimits = xlim;
end
if ~isempty(ylimits)
    yl = ylim;
    ylimits(ylimits == -inf) = yl(1);
    ylimits(ylimits == inf) = yl(2);
    ylim(ylimits);
else
    ylimits = ylim;
end

switch logScale
    case 0
        xrange = xlimits(2)-xlimits(1);
        yrange = ylimits(2)-ylimits(1);
    case 1
        xrange = log(xlimits(2))-log(xlimits(1));
        yrange = ylimits(2)-ylimits(1);
    case 2
        xrange = xlimits(2)-xlimits(1);
        yrange = log(ylimits(2))-log(ylimits(1));
    case 3
        xrange = log(xlimits(2))-log(xlimits(1));
        yrange = log(ylimits(2))-log(ylimits(1));
end

for i = 1:nLines
    %candidateTextLocations{i} = markerLocations{i};
    candidateTextLocations{i} = zeros(1,0);
    for j = 1:length(x{i}) % Restrict to values on line that are contained in the plot
        if x{i}(j) >= xlimits(1) && x{i}(j) <= xlimits(2) && ...
                y{i}(j) >= ylimits(1) && y{i}(j) <= ylimits(2)
            candidateTextLocations{i}(1,end+1) = j;
        end
    end
end

% Now either label the lines, or just make a legend
if ~isempty(legendStr)
    if labelLines % Try to label the lines
        for i = 1:nLines
            
            switch labelLoc
                case 0 % Random Placement on the line
                    m = candidateTextLocations{i}(ceil(rand*length(candidateTextLocations{i})));
                case 1 % As far as possible from nearest currently-placed label
                    if i == 1
                        m = candidateTextLocations{i}(ceil(rand*length(candidateTextLocations{i})));
                    else
                        furthest = 0;
                        for m = candidateTextLocations{i}
                            d = figDist(x{i}(m),y{i}(m),textLocations(1:i-1,1),textLocations(1:i-1,2));
                            closest = min(d);
                            if closest > furthest
                                furthest = closest;
                                bestM = m;
                            end
                        end
                        %fprintf('Furthest from %d is %d\n',i,bestM);
                        m = bestM;
                    end
                case 2 % Place label as far as possible from other lines and figure boundary
                    furthest = 0;
                    bestM = candidateTextLocations{i}(ceil(rand*length(candidateTextLocations{i})));
                    for m = candidateTextLocations{i}
                        closest = inf;
                        
                        % Compare to figure boundaries
                        d1 = figDist(x{i}(m),y{i}(m),x{i}(m),ylimits(1));
                        d2 = figDist(x{i}(m),y{i}(m),x{i}(m),ylimits(2));
                        d3 = figDist(x{i}(m),y{i}(m),xlimits(1),y{i}(m));
                        d4 = figDist(x{i}(m),y{i}(m),xlimits(2),y{i}(m));
                        closest = min([d1 d2 d3 d4]);
                        
                        % Compare to other lines
                        for j = [1:i-1 i+1:nLines]
                            d = figDist(x{i}(m),y{i}(m),x{j},y{j});
                            closest = min([closest d]);
                        end
                                                
                        if closest > furthest
                            furthest = closest;
                            bestM = m;
                        end
                    end
                    m = bestM;
            end
            textLocations(i,:) = [x{i}(m),y{i}(m)];
            
            % Draw text
            stringVal = [repmat(' ',1,markerSize/4),legendStr{i}];
            stringVal = legendStr{i};
            hText(i) = text(x{i}(m),y{i}(m),...
                stringVal,...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','center',...
                'FontSize',12,'FontWeight','bold');
            
            % Fill background and add edge
            set(hText(i),'BackgroundColor',[1 1 1]);
            set(hText(i),'EdgeColor',[0 0 0]);
            
            % Add Color
            if ~isempty(colors)
                if iscell(colors)
                    set(hText(i),'Color',colors{1+mod(i-1,length(colors))});
                    %set(hText(i),'EdgeColor',colors{1+mod(i-1,length(colors))});
                else
                    set(hText(i),'Color',colors(1+mod(i-1,size(colors,1)),:));
                   % set(hText(i),'EdgeColor',colors(1+mod(i-1,size(colors,1)),:));
                end
            end
            
            % Change Edge Style
            if ~isempty(lineStyles) && useLines
                if ~isempty(lineStyles)
                    if ~isempty(lineStyles{1+mod(i-1,length(lineStyles))})
                        %set(hText(i),'LineStyle',lineStyles{1+mod(i-1,length(lineStyles))});
                    end
                end
            end
            
            if labelRotate
                p1 = max(m-1,1);
                p2 = min(m+1,length(x{i}));
                slope = figDistY(y{i}(p2),y{i}(p1))/figDistX(x{i}(p2),x{i}(p1));
                set(hText(i),'Rotation',180*atan(slope)/pi);
            end
        end
    else % Standard legend
        h = legend(h,legendStr);
        set(h,'FontSize',10,'FontWeight','normal');
        set(h,'Location','Best');
        set(h,'Location',legendLoc);
    end
end

set(gcf, 'PaperPositionMode', 'auto');


% Changes the width/style/marker of a line
    function [h] = applyStyle(h,i,lineWidth,colors,lineStyles,markers,markerSpacing)
        hold on;
        set(h,'LineWidth',lineWidth);
        if ~isempty(colors)
            if iscell(colors)
                set(h,'Color',colors{1+mod(i-1,length(colors))});
            else
                set(h,'Color',colors(1+mod(i-1,size(colors,1)),:));
            end
        end
        if ~isempty(lineStyles) && useLines
            if isempty(lineStyles)
                set(h,'LineStyle','-');
            else
                if ~isempty(lineStyles{1+mod(i-1,length(lineStyles))})
                    set(h,'LineStyle',lineStyles{1+mod(i-1,length(lineStyles))});
                end
            end
        end
        if ~isempty(markers)
            if ~isempty(markers{1+mod(i-1,length(markers))})
                if isempty(markerSpacing)
                    set(h,'Marker',markers{1+mod(i-1,length(markers))});
                    set(h,'MarkerSize',markerSize);
                    if fillFace
                        set(h,'MarkerFaceColor',[1 1 .9]);
                    end
                end
            end
        end
    end

% Computes distance on figure, taking into account range of x/y-values
% and use of log-scale
    function [d] = figDist(x1,y1,x2,y2)
        v = [figDistX(x1,x2(:)')
            figDistY(y1,y2(:)');];
        d = sqrt(sum(v.^2));
    end

    function [dX] = figDistX(x1,x2)
        xlimits = xlim;
        switch logScale
            case {0,2}
                dX = (x1-x2)./(xlimits(2)-xlimits(1));
            case {1,3}
                dX = (log(x1)-log(x2))./(log(xlimits(2))-log(xlimits(1)));
        end
    end
    function [dY] = figDistY(y1,y2)
        ylimits = ylim;
        switch logScale
            case {0,1}
                dY = (y1-y2)/(ylimits(2)-ylimits(1));
            case {2,3}
                dY = (log(y1)-log(y2))/(log(ylimits(2))-log(ylimits(1)));
        end
    end
end