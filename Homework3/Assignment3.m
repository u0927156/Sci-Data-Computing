%% Read in Filip data set
clear; clc;
fileID = fopen('NIST_Filip_trim.txt', 'r');

Input = fscanf(fileID, '%f %f', [2 inf])';

fclose(fileID);

%% Sort input
Input = sortrows(Input,2);


x = Input(:,2);
y = Input(:,1);

figure(301)
plot(x,y, 'x')
xlabel('x'); ylabel('y'); title('NIST Filip Dataset')
saveas(gcf, 'NISTDataset.png')

%% 
full_vander = vander(x);

ms = [9 11 13 15 17 19 21 23 25];


normal_coefficients_cell = cell(1,length(ms));
qr_coefficients_cell = cell(1,length(ms));

for n = 1:length(ms)
    m = ms(n);
    
    A = full_vander(:, (end-m):end);

    yy = A'*y;
    B = A'*A;
    normal_coefficients = B\yy; % matlab gives ill conditioned warning
    normal_coefficients_cell{n} = normal_coefficients;
    
    norm(A*normal_coefficients-y,2);

    [Q, R] = qr(A);

    yy = Q'*y;
    W = Q'*Q;
    qr_coefficients = R\yy;
    qr_coefficients_cell{n} = qr_coefficients;
    norm(A*qr_coefficients-y,2);
end 


figure(302)
plot(x,y, 'x')
xlabel('x'); ylabel('y'); title('NIST Filip Dataset and Solutions')
to_plot = [1 4  9];


x_plot = linspace(min(x), max(x), 1001);
vander_plot = vander(x_plot);
hold on
legend_info = cell(1, 2*length(to_plot)+1);
legend_info{1} = 'Normal Data';
cell_counter = 2;
for n = to_plot
    m = ms(n);
    A = vander_plot(:, (end-m):end);
    
    
    normal_coefficients = normal_coefficients_cell{n};
    normal_to_plot = A*normal_coefficients;
    legend_info{cell_counter} = ['Normal ' num2str(m) ' Coefficients'];
    cell_counter = cell_counter+1;
    
    qr_coefficients = qr_coefficients_cell{n};
    qr_to_plot = A*qr_coefficients;
    legend_info{cell_counter} = ['QR ' num2str(m) ' Coefficients'];
    cell_counter = cell_counter+1;
    
    plot(x_plot, normal_to_plot)
    
    plot(x_plot, qr_to_plot)
end
ylim([0.76 0.94])
hold off
legend(legend_info, 'location', 'southeast')
saveas(gcf, 'NISTDatasetSolutions.png')
%% 1c Calculate least square errors 
full_vander = vander(x);

normal_error = zeros(1, length(ms));
qr_error = zeros(1, length(ms));
for n = 1:length(ms)
    m = ms(n);
    A = full_vander(:, (end-m):end);
    
    
    normal_coefficients = normal_coefficients_cell{n};
    normal_for_error = A*normal_coefficients;
    normal_error(n) = sum((y-normal_for_error).^2);
    
    qr_coefficients = qr_coefficients_cell{n};
    qr_for_error = A*qr_coefficients;
    qr_error(n) = sum((y-qr_for_error).^2);
end

figure(303)
semilogy(ms, normal_error,'--x')

hold on
semilogy(ms, qr_error, '-.o')
hold off

title('Least Squares Error for NIST Dataset')
xlabel('Polynomial Coefficients'); ylabel('Least Squares Error')
xticks(ms)
legend('Normal Equations', 'QR Method', 'location', 'northwest')
saveas(gcf, 'NISTDatasetError.png')
%% Scaling
y_scaled = 2/(max(y)-min(y)).*y - (max(y)+min(y))/(max(y)-min(y));
x_scaled = 2/(max(x)-min(x)).*x - (max(x)+min(x))/(max(x)-min(x));

figure(304)
plot(x_scaled, y_scaled, 'x')
title('Rescaled Data')

full_vander = vander(x_scaled);

ms = [9 11 13 15 17 19 21 23 25];


normal_coefficients_cell = cell(1,length(ms));
qr_coefficients_cell = cell(1,length(ms));

for n = 1:length(ms)
    m = ms(n);
    
    A = full_vander(:, (end-m):end);

    yy = A'*y_scaled;
    B = A'*A;
    normal_coefficients = B\yy; % matlab gives ill conditioned warning
    normal_coefficients_cell{n} = normal_coefficients;
    
    norm(A*normal_coefficients-y_scaled,2);

    [Q, R] = qr(A);

    yy = Q'*y_scaled;
    W = Q'*Q;
    qr_coefficients = R\yy;
    qr_coefficients_cell{n} = qr_coefficients;
    norm(A*qr_coefficients-y_scaled,2);
end 


figure(305)
plot(x_scaled,y_scaled, 'x')
xlabel('x'); ylabel('y'); title('Rescaled NIST Filip Dataset')
saveas(gcf, 'ScaledNISTDataset.png')

to_plot = [1 4 9];

x_plot = linspace(min(x_scaled), max(x_scaled), 1001);
vander_plot = vander(x_plot);
hold on
legend_info = cell(1, 2*length(to_plot)+1);
legend_info{1} = 'Normal Data';
cell_counter = 2;


for n = to_plot
    m = ms(n);
    A = vander_plot(:, (end-m):end);
    
    
    normal_coefficients = normal_coefficients_cell{n};
    normal_to_plot = A*normal_coefficients;
    legend_info{cell_counter} = ['Normal ' num2str(m) ' Coefficients'];
    cell_counter = cell_counter+1;
    
    qr_coefficients = qr_coefficients_cell{n};
    qr_to_plot = A*qr_coefficients;
    legend_info{cell_counter} = ['QR ' num2str(m) ' Coefficients'];
    cell_counter = cell_counter+1;
    
    plot(x_plot, normal_to_plot)
    
    plot(x_plot, qr_to_plot)
end
ylim([-1 1])
hold off
legend(legend_info, 'location', 'southeast')
saveas(gcf, 'ScaledNISTDatasetSolution.png')
%% Calculate least square errors for scaled data
full_vander = vander(x_scaled);

normal_error = zeros(1, length(ms));
qr_error = zeros(1, length(ms));
for n = 1:length(ms)
    m = ms(n);
    A = full_vander(:, (end-m):end);
    
    
    normal_coefficients = normal_coefficients_cell{n};
    normal_for_error = A*normal_coefficients;
    normal_error(n) = sum((y_scaled-normal_for_error).^2);
    
    qr_coefficients = qr_coefficients_cell{n};
    qr_for_error = A*qr_coefficients;
    qr_error(n) = sum((y_scaled-qr_for_error).^2);
end

figure(306)
semilogy(ms, normal_error,'--x')

hold on
semilogy(ms, qr_error, '-.o')
hold off

title('Least Squares Error for NIST Dataset')
xlabel('Polynomial Coefficients'); ylabel('Least Squares Error')
xticks(ms)
legend('Normal Equations', 'QR Method', 'location', 'northwest')
saveas(gcf, 'ScaledNISTDatasetError.png')
%% SVD Image Processing
clear; clc

Gatlin = 1;

if Gatlin
    load gatlin
    filename_base = 'gatlin';
else
    load durer
    filename_base = 'durer';
end

figure(307)
image(X), colormap('gray')
title('Full Image','fontsize',16)
axis off
saveas(gcf, [filename_base '_Base.png'])
    

rs = [2 8 32 128];

svd_time = zeros(1,length(rs));
rsvd_time = zeros(1,length(rs));

figure_counter = 310;
for n = 1:length(rs)
    r = rs(n);
    
    
    tic 
    [U,S,V] = svd(X);
    Xk = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    svd_time(n) = toc;
    
    
    tic
    [U_r,S_r,V_r] = rsvd(X,r);
    Xr = U_r*S_r*V_r';
    rsvd_time(n) = toc;
    
    rsvd_singular_values{n} = diag(S_r);
    
    figure(figure_counter)
    image(Xk), colormap(map)
    title(sprintf('Rank-%d approximation with SVD',r),'fontsize',16)
    axis off
    saveas(gcf, [filename_base '_SVD' num2str(r) '.png'])
    
    figure(figure_counter+1)
    image(Xr), colormap(map)
    title(sprintf('Rank-%d approximation with rSVD',r),'fontsize',16)
    axis off
    saveas(gcf, [filename_base '_rSVD' num2str(r) '.png'])
    
    
    figure_counter =figure_counter +2;
end

figure(figure_counter)

semilogy(diag(S), 'LineWidth', 2)
title('Singular Values of SVD and rSVD')
ylabel('Magnitude of Singular Value')
xlabel('Value')
hold on
for n = 1:length(rs)
    semilogy(rsvd_singular_values{n},'x')
end

xlim([0, 128])
legend('SVD', 'rSVD, r=2','rSVD, r=8','rSVD, r=32','rSVD, r=128',...
    'location', 'northeast')
hold off
saveas(gcf, [filename_base '_SingularValueComparison.png'])


%% rSVD and Power Method
clear; clc;

Gatlin = 0;

if Gatlin
    load gatlin
    filename_base = 'gatlin';
else
    load durer
    filename_base = 'durer';
end
r = 32;
qs = [1 2 3];

rsvd_time = zeros(1,length(qs));

figure_counter = 330;
for n = 1:length(qs)
    q = qs(n);

    % calculate time to find rSVD and reconstruct matrix.
    tic
    [U_r,S_r,V_r] = power_rsvd(X,r, q); % use modified version of rSVD
    Xr = U_r*S_r*V_r';
    rsvd_time(n) = toc;
    
    power_rsvd_singular_values{n} = diag(S_r);

    
    figure(figure_counter)
    image(Xr), colormap(map)
    title(sprintf('Rank-%d rSVD, q = %d',r, q),'fontsize',16)
    axis off
    saveas(gcf, [filename_base '_power_rSVD' num2str(q) '.png'])
    
    
    figure_counter =figure_counter +1;
end

figure(figure_counter)
[U,S,V] = svd(X);
hold off
relevent_s = diag(S(1:r, 1:r));
semilogy(diag(S(1:r, 1:r)))
hold on

error_norms = zeros(1,length(qs));
for n = 1:length(qs)
    semilogy(power_rsvd_singular_values{n},'x')
    error_norms(n) = norm(relevent_s-power_rsvd_singular_values{n},2);
end
hold off
legend('SVD', 'q=1', 'q=2', 'q=3')

error_norms

%% Matrix
clear; clc
X = [611.0 196.0 -192.0 407.0 -8.0 -52.0 -49.0 29.0;
196.0 899.0 113.0 -192.0 -71.0 -43.0 -8.0 -44.0;
192.0 113.0 899.0 196.0 61.0 49.0 8.0 52.0;
407.0 -192.0 196.0 611.0 8.0 44.0 59.0 -23.0;
-8.0 -71.0 61.0 8.0 411.0 -599.0 208.0 208.0;
-52.0 -43.0 49.0 44.0 -599.0 411.0 208.0 208.0;
-49.0 -8.0 8.0 59.0 208.0 208.0 99.0 -911.0;
29.0 -44.0 52.0 -23.0 208.0 208.0 -911.0 99.0];

errors_svd = zeros(1, 8);
for r = 1:8
    [U,S,V] = svd(X);
    Xk = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    errors_svd(r)= norm(Xk-X,'fro');
end

errors_svd 


errors_rsvd = zeros(1, 8);
for r = 1:8
    [U_r,S_r,V_r] = rsvd(X,r);
    Xr = U_r*S_r*V_r';
    errors_rsvd(r)= norm(Xr-X,'fro');
end
errors_rsvd

errors_rsvd_power = zeros(1, 8);
for r = 1:8
    [U_r,S_r,V_r] = power_rsvd(X,r,3);
    Xr = U_r*S_r*V_r';
    errors_rsvd_power(r)= norm(Xr-X,'fro');
end
errors_rsvd_power


figure(400)
semilogy(errors_svd, '-x')
hold on 
semilogy(errors_rsvd, '--o')
semilogy(errors_rsvd_power, '-.*')
hold off

title('Error of SVD and rSVD')
xlabel('Number of Singular Values')
ylabel('Error, Frobenius Norm')
legend('SVD', 'rSVD', 'rSVD, q=3',...
    'location', 'southwest')

saveas(gcf, 'TrickyMatrix.png')


