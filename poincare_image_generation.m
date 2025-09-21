%diff poincare image
xpdx = RR_Interval(1:end-1);
xmdy = RR_Interval(2:end);
figure
plot(xpdx, xmdy, 'b', 'LineWidth',1.5 )
xlim([-0.8 0.8])
ylim([-0.8 0.8])
frame = getframe;
image_matrix = frame2im(frame); % Convert the frame to an image matrix
resized_image = imresize(image_matrix, [224, 224]);
image_features; % function for image domain feature extraction