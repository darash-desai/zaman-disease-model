function ploterrors(X, Y, color)
    %% Plots a data set with error region of +/- 2 std.   
    %
    % If a vector is provided for the Y data, the data is simply plotted
    % without any error region.
    %
    % Parameters:
    %   - X: The X dataset.
    %   - Y: Matrix of Y data where each row represents replicate data.   
      
    if (size(Y, 1) == 1)
        plot(X, Y, 'Color', color, 'LineWidth', 2);
        return;
    end

    avg = mean(Y, 1);
    error = 2 * std(Y);
    pos_error = avg + error;
    neg_error = avg - error;   

    % Turn hold on if it is not already on.
    isHoldOn = ishold;
    if ~isHoldOn
        hold on;
    end    

    X2 = [X, fliplr(X)];
    area = [pos_error, fliplr(neg_error)];
    fill(X2, area, 'r', 'FaceAlpha', 0.5, 'FaceColor', color, 'LineStyle', 'none', 'HandleVisibility', 'off');    
    plot(X, avg, 'Color', color, 'LineWidth', 2);
    ylim([0 inf]);

    % Return hold to it's previous off state.
    if ~isHoldOn
        hold off;
    end
end