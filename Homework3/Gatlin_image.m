close all

% Select to load gatlin or durer image with this variable
Gatlin = 0;

if Gatlin
    load gatlin
    filename_base = 'gatlin';
else
    load durer
    filename_base = 'durer';
end




% load the "gatlin" image data, built-in to MATLAB

[U,S,V] = svd(X);
 
% "gatlin" stores the image as the variable "X"

figure(1),clf
 % plot the singular values

semilogy(diag(S),'b.','markersize',20)


set(gca,'fontsize',16)

title('singular values of the "gatlin" image matrix')

xlabel('k'), ylabel('\sigma_k')

figure(2),clf

% plot the original image

image(X), colormap(map)

% image: MATLAB command to display a matrix as image

axis equal, axis off
 
title('true image (rank 480)','fontsize',16)



 
% plot the optimal rank-k approximation for each relevant k value
new_figure = 3;
for k = [2 8 32 128]
    
    figure(new_figure),clf
    new_figure = new_figure+1;
    Xk = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';

    image(Xk), colormap(map)

    axis equal, axis off

    title(sprintf('best rank-%d approximation',k),'fontsize',16)

    saveas(gcf, [filename_base num2str(k) '.png'])
end


