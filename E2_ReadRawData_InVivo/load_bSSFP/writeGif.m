%% write GIF

h = figure('Color','White','Position',[300 300 600 500]);
axis tight manual equal % this ensures that getframe() returns a consistent size
filename = 'G.gif';
for n = 1:size(Gimg,3)
    % Draw plot for y = x.^n
    imagesc(Gimg(:,:,n),[0 3]); axis equal;axis off;  %colorbar; % 12th slice
    drawnow
      % Capture the plot as an image
      frame = getframe(h);
      im = frame2im(frame);
     
      [imind,cm] = rgb2ind(im,256);
      % Write to the GIF File
      if n == 1
          imwrite(imind,cm,filename,'gif', 'DelayTime', 0.1, 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif', 'DelayTime', 0.1,'WriteMode','append');
      end
end
close(h);