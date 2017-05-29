function orbitPlotter(x, y)

    colors = [204 0 153; 255 0 0; 0 0 0]/255;
    colorsRepeated = ones(length(x), 3);
    sizes = [20; 5; 2];
    sizesRepeated = ones(length(x), 1);

    for i=1:length(x)
        sizesRepeated(i) = sizes(mod(i-1, 3)+1);

        for j=1:3
            colorsRepeated(i,j) = colors(mod(i-1,3)+1,j);
        end

    end

    scatter(x, y, sizesRepeated, colorsRepeated, 'filled')

    xlabel('x');
    ylabel('y');

    hold on;
    h = zeros(3, 1);
    for i = 1:3
        colors(i,:)
        h(i) = scatter(NaN, NaN, 20, colors(i,:), 'filled');
    end
    legend(h, 'Sun','Jupiter', 'lol');
    hold off;
end

