function orbitPlotter(x, y, rotate)

    copyX = x;
    copyY = y;

    colors = [255 0 0; 255, 119, 51; 0 0 0]/255;
    sizes = [40; 20; 100];
    markers = ['*', 'o', '.'];
    nBodies = 3;

    hold on;
    
%     if rotate
%         for i=1:length(x)/nBodies
%             rotAngle = -atan2(y((i-1)*nBodies+2), x((i-1)*nBodies+2));
%             for j=1:nBodies
%                 copyX(i) = cos(rotAngle)*x((i-1)*nBodies+j) - sin(rotAngle)*y((i-1)*nBodies+j);
%                 copyY(i) = sin(rotAngle)*x((i-1)*nBodies+j) + cos(rotAngle)*y((i-1)*nBodies+j);
%             end
%         end
%     end
        
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
    end
    
    legend('Sun','Jupiter', 'test mass');
    
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
    pbaspect([1 1 1]);
    hold off;
end