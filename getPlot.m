function[] = getPlot(x, y, legended)
    N = size(y);

    figure
    hold on
    grid on
    for i=1:N(1)
        plot(x(:), y(i, :), 'DisplayName', strcat(legended, "_", num2str(i)));
    end
    xlabel('time, s');
    ylabel(legended);
    hold off
    legend
end