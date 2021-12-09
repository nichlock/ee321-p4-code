clear;clc;close all;
        
iterations = 40;

N = 100000; % Number of symbols 

plot_n = [29, 32];

max_bound = 2.5;


dist = @(cpx1, cpx2) (sqrt((real(cpx2 - cpx1) .^ 2) + ((imag(cpx2 - cpx1)) .^ 2)));

alfs_ = 14; % Axes Label Font Size

s2 = sqrt(2);
schemes = [
    1, 0;
    0, 1;
    -1, 0;
    0, -1;
    1./s2, 1./s2;
    -1./s2, 1./s2;
    -1./s2, -1./s2;
    1./s2, -1./s2;
    1+1i, 1;
    ];

for state = 1:9

    x_r_fig = figure('Position', [0, 500-(state*50), 600, 500]);
    x_r_err_fig = figure('Position', [600, 500-(state*50), 600, 500]);
    const_fig = figure('Position', [1200, 500-(state*50), 600, 500]);
%     sgtitle(['$x$ and $r$ for Scheme $s_m = (', num2str(schemes(state, 1)), ', ', num2str(schemes(state, 2)), ')$'],...
%         'interpreter','latex');
%     ax_const = subplot(3, 4, [3:4, 7:8, 11:12]);
%     ax_re_x_r = subplot(3, 4, [1, 5]);
%     ax_im_x_r = subplot(3, 4, [2, 6]);

    figure(x_r_err_fig);
    ax_re_x_r = subplot(3, 2, [1, 3]);
    ax_im_x_r = subplot(3, 2, [2, 4]);
    ax_err = subplot(3, 2, [5, 6]);

    err_rates = zeros(iterations + 1, 1);
    snrs = zeros(iterations + 1, 1);

    for it = 0:iterations
        snr = it - (iterations ./ 2);
        snrs(it+1) = snr;

        var_n = 0.01; % Noise variance
        A = sqrt(real((exp(snr)) ./ 10) .* (var_n^2));
        x = A * ((2 * (randn(N,1)>0)) - 1);
        x_im = A * ((2 * (randn(N,1)>0)) - 1) * 1i;
        x = (x .* schemes(state, 1)) + (x_im .* schemes(state, 2));

        x_possibilities = unique(x);

        x_known = zeros(N, 1);
        for x_ind = 1:length(x)
            for guess_ind = 1:length(x_possibilities)
                if(x(x_ind) == x_possibilities(guess_ind))
                    x_known(x_ind) = guess_ind;
                end
            end
        end


        n = randn(N,1)*sqrt(var_n); % Zero-mean Gaussian noise with variance var_n
        n_im = randn(N,1)*sqrt(var_n) * 1i; % Zero-mean Gaussian noise with variance var_n
        r = x + n + n_im; % Received signal

        % Using bounds, finds error of r to x
        r_arr = zeros(1, 2);
        r_guesses = zeros(N, 1);
        for r_ind = 1:length(r)
            top_guess_ind = 0;
            r_guess_val = +inf;
            for guess_ind = 1:length(x_possibilities)
                guess_distance = dist(r(r_ind), x_possibilities(guess_ind));
                if(r_guess_val > guess_distance)
                    top_guess_ind = guess_ind;
                    r_guess_val = guess_distance;
                end
            end
            r_guesses(r_ind) = top_guess_ind;
        end
        success_rate = r_guesses == x_known;

        error_rates = 1 - success_rate;
        err_rates(it + 1) = sum(error_rates) ./ length(error_rates);
        
        if(ismember(it, plot_n))
            figure(x_r_fig);
            n_to_plot = 1000;
            sgtitle('Real and Imaginary Components of $x$ and $r$', 'interpreter', 'latex');
            subplot(2, 2, 1);
            hold on;
            scatter(1:n_to_plot, real(x(1:n_to_plot)));
            title('$Re[x]$', 'interpreter', 'latex');
            subplot(2, 2, 2);
            hold on;
            scatter(1:n_to_plot, real(r(1:n_to_plot)));
            title('$Re[r]$', 'interpreter', 'latex');

            subplot(2, 2, 3);
            hold on;
            scatter(1:n_to_plot, imag(x(1:n_to_plot)));
            title('$Im[x]$', 'interpreter', 'latex');
            subplot(2, 2, 4);
            hold on;
            scatter(1:n_to_plot, imag(r(1:n_to_plot)));
            title('$Im[r]$', 'interpreter', 'latex');
        end
    
        if(ismember(it, plot_n))
            figure(x_r_err_fig);
            subplot(ax_re_x_r);
            hold(ax_re_x_r, 'on');
            scatter(1:N, real(r));
            scatter(1:N, real(x));
            axis([1, N, -max_bound, max_bound]);
            title('Real $x$ and $r$', 'interpreter','latex', 'FontSize',alfs_);
            xlabel('$n$', 'interpreter','latex', 'FontSize',alfs_);
            ylabel('$Re[r, x]$', 'interpreter','latex', 'FontSize',alfs_);
            
            subplot(ax_im_x_r);
            hold(ax_im_x_r, 'on');
            scatter(1:N, imag(r));
            scatter(1:N, imag(x));
            axis([1, N, -max_bound, max_bound]);
            title('Imaginary $x$ and $r$', 'interpreter','latex', 'FontSize',alfs_);
            xlabel('$n$', 'interpreter','latex', 'FontSize',alfs_);
            ylabel('$Im[r, x]$', 'interpreter','latex', 'FontSize',alfs_);
            
            figure(const_fig);
            %subplot(ax_const);
            hold(gca, 'on')
            scatter(real(r), imag(r), 50, '.');
            scatter(real(x), imag(x), 100, '.', 'MarkerFaceColor', '#D95319');
            if(it == plot_n(end))
                xlabel('$Re[r, x]$', 'interpreter','latex', 'FontSize',alfs_);
                ylabel('$Im[r, x]$', 'interpreter','latex', 'FontSize',alfs_);
            
                snr_str = strjoin(cellstr(num2str(snrs(plot_n(:) + 1))),',');
                title(['Constellation for $x$ and $r$, $SNR=', snr_str, '$ dB'], 'interpreter','latex', 'FontSize',alfs_);
                
                for x_ind = 1:(length(x_possibilities))
                    for x_ind2 = 1:(length(x_possibilities))
                        if((x_ind == x_ind2))
                            continue;
                        end
                        bound = real(x_possibilities(x_ind2) - x_possibilities(x_ind))...
                               / imag(x_possibilities(x_ind2) - x_possibilities(x_ind));
                        if(bound == inf)
                            bound = 1000000;
                        end
                        if(bound == 1 || bound == -1)
                            continue
                        end
                        plot(-10:10, (-10:10) .* -bound, '--', 'LineWidth', 2, 'Color', '#77AC30');
                    end
                end
                
                axis([-max_bound, max_bound, -max_bound, max_bound]);
                
                [~, leg_obj] = legend(...
                    ['$r$, $SNR=', num2str(snrs(plot_n(1) + 1)), '$dB'],...
                    '$x$',...
                    ['$r$, $SNR=', num2str(snrs(plot_n(2) + 1)), '$dB'],...
                    '$x$',...
                    'Boundries',...
                    'FontSize',10, 'Interpreter','latex');
                set(findobj(leg_obj, 'type', 'patch'), 'Markersize', 20);
                if(state > 4)
                    set(findobj(leg_obj, '-property', 'Location'), 'Location', 'north');
                end
                drawnow update
            end

        end
        if(err_rates(it+1) == 0 && it >= plot_n(end))
            break;
        end
    end

    figure(x_r_err_fig);
    ax = subplot(ax_err);
    semilogy(snrs, err_rates);
    title('Error Rate vs. SNR', 'interpreter','latex', 'FontSize',alfs_);
    axis([-20, 20, 0, 1]);
    drawnow update
    ax.Position = [0.13,0.11,0.775,0.17];
    drawnow update

    
    %saveas(x_r_fig,['xr/s_', num2str(state), '.png'], 'png')
    %saveas(x_r_err_fig,['err/s_', num2str(state), '.png'], 'png')
    %saveas(const_fig,['constellation/s_', num2str(state), '.png'], 'png')

    err_thresh = 0.01;
    for it = 1:iterations
        if((err_rates(it) <= err_thresh) || it == iterations)
            if(it == iterations)
                fprintf('Scheme = %d; error < %2.2f%% not found\n', it, err_thresh*100);
                break;
            end
            fprintf('Scheme = %d; Error = %2.2f%% of signal at SNRs = %d dB\n', state, err_rates(it)*100, snrs(it));
            break;
        end
    end
end
