function orbitPlotter(x, y)

    colors = [255 0 0; 255, 119, 51; 0 0 0]/255;
    sizes = [40; 20; 2];
    markers = ['*', 'o', '.'];
    nBodies = 3;

    hold on;
    
    for i=1:nBodies
        bodyX = x(i:nBodies:length(x));
        bodyY = y(i:nBodies:length(y));
        scatter(bodyX, bodyY, sizes(i), colors(i,:), markers(i))
    end
    
    legend('Sun','Jupiter', 'test mass');
    
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
    pbaspect([1 1 1]);
    hold off;
end