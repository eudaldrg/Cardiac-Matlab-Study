function varargout = plotshaded(x,y,fstr);
 
if size(y,1)>size(y,2)
    y=y';
end;

% if size(y,1)==1 % just plot one line
%     plot(x,y,'fstr');
% end;
 
if size(y,1)==2 %plot shaded area
    px=[x,fliplr(x)]; % make closed patch
    py=[y(1,:), fliplr(y(2,:))];
    patch(px,py,1,'FaceColor',fstr,'EdgeColor','none');
end;
 
% if size(y,1)==3 % also draw mean
%     px=[x,fliplr(x)];
%     py=[y(1,:), fliplr(y(3,:))];
%     patch(px,py,1,'FaceColor',fstr,'EdgeColor','none');
%     plot(x,y(2,:),fstr);
% end;
 
alpha(.2); % make patch transparent

