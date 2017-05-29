function orbitPlotter(x, y, rotate, limits)

    copyX = x;
    copyY = y;

    colors = [255 0 0; 255, 119, 51; 0 0 0]/255;
    sizes = [40; 20; 4];
    markers = ['*', 'o', '.'];
    nBodies = 3;

    for i=1:nBodies
        bodyX = copyX(i:nBodies:length(x));
        bodyY = copyY(i:nBodies:length(y));
        
        if rotate
            if i==1
                rotAngles = -atan2(bodyY, bodyX)-pi;
            end
            
            for j=1:length(bodyX)
                newX = cos(rotAngles(j))*bodyX(j) - sin(rotAngles(j))*bodyY(j);
                newY = sin(rotAngles(j))*bodyX(j) + cos(rotAngles(j))*bodyY(j);
                bodyX(j) = newX;
                bodyY(j) = newY;
            end
        end
        
        scatter(bodyX, bodyY, sizes(i), colors(i,:), markers(i))
        hold on;
    end
    
    legend('Sun','Jupiter', 'test mass');
    
    if limits
        xlim([-1.1 1.1]);
        ylim([-1.1 1.1]);
    end
    
    pbaspect([1 1 1]);
    hold off;
end