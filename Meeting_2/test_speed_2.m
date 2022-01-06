clear
% close all
clc

addpath([pwd '\Integrator']);

%% settings

noise_factors = [1:7];
use_feedback = 0;

n = 1; % number of trials

res_name = ['compare4_noise_factor_' num2str(noise_factors(1)) '_',...
    num2str(noise_factors(end)) '_' 'trials_' num2str(n)];


%% initialise result storage
time_3 = nan(n,1);
time_2 = nan(n,1);

Results3 = cell(n,1);
Results2 = cell(n,1);


%% Run both versions of the code for the given noise factors...
% ... and repeat n times

log_folder = ['logfiles_' res_name];
mkdir(log_folder);

for i=1:n

    diary([pwd '/' log_folder '/' 'example2_2c_i' num2str(i) '.txt'])
    t0 = tic;
    [R1] = example2_2c(0,noise_factors,0.167,use_feedback);
    t1 = toc(t0);
    Results3{i} = R1;
    time_3(i) = t1;
    diary off


    diary([pwd '/' log_folder '/' 'example2_2b_i' num2str(i) '.txt'])
    t0 = tic;
    [R2] = example2_2b(0,noise_factors,0.167,use_feedback);
    t2 = toc(t0);
    Results2{i} = R2;
    time_2(i) = t2;
    diary off

end

%% save results

save([res_name '.mat'],'time_3','time_2','Results3','Results2');

%% load results

load([res_name '.mat'])

%% plot comparison

% noise_factors = [1:6];

cs1 = linspecer(length(noise_factors)+1,'blue');
cs1 = cs1(2:end,:);

cs2 = linspecer(length(noise_factors)+1,'green');
cs2 = cs2(2:end,:);

cs3 = linspecer(length(noise_factors)+1,'red');
cs3 = cs3(2:end,:);

csd = linspecer(length(noise_factors)+1,'grey');
csd = csd(2:end,:);

nv = 2;
nh = 6;

scs = get(0,'ScreenSize');

figure('Position',[0,40,scs(3),scs(4)-120]);

subplot(nv,nh,1)
pl1=plot(1:n,time_3,'o','Color',cs3(end,:),'MarkerFaceColor',cs3(end,:),'DisplayName','elapsed time');
hold on
pl2=plot(1:n,time_2,'d','Color',cs2(end,:),'MarkerFaceColor',cs2(end,:),'DisplayName','elapsed time');
xlabel('trial')
ylabel('time (s)')
title('Total trial duration')
axis tight;
yl = get(gca, 'ylim');
ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]);
xlim([0.5,n+0.5]);

for i=1:n
    R3i = Results3{i};
    R2i = Results2{i};

    t1_tot_i =0;
    t2_tot_i =0;

    for j=1:length(noise_factors)
        R3ij = R3i.(['R' num2str(noise_factors(j))]);
        R2ij = R2i.(['R' num2str(noise_factors(j))]);

        subplot(nv,nh,3)
        plot(i,R3ij.stats.t_wall_total,'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
        hold on
        plot(i,R2ij.stats.t_wall_total,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('trial')
        ylabel('time (s)')
        title('Total CPU time')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); xlim([0.5,n+0.5]);
        t1_tot(i,j) = R3ij.stats.t_wall_total;
        t2_tot(i,j) = R2ij.stats.t_wall_total;

        t1_tot_i = t1_tot_i + R3ij.stats.t_wall_total;
        t2_tot_i = t2_tot_i + R2ij.stats.t_wall_total;


        subplot(nv,nh,4)
        t_feval_1ij = R3ij.stats.t_wall_nlp_f + R3ij.stats.t_wall_nlp_g + ...
            R3ij.stats.t_wall_nlp_grad + R3ij.stats.t_wall_nlp_grad_f + ...
            R3ij.stats.t_wall_nlp_hess_l + R3ij.stats.t_wall_nlp_jac_g;
        t_feval_2ij = R2ij.stats.t_wall_nlp_f + R2ij.stats.t_wall_nlp_g + ...
            R2ij.stats.t_wall_nlp_grad + R2ij.stats.t_wall_nlp_grad_f + ...
            R2ij.stats.t_wall_nlp_hess_l + R2ij.stats.t_wall_nlp_jac_g;
        plot(i,t_feval_1ij,'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
        hold on
        plot(i,t_feval_2ij,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('trial')
        ylabel('time (s)')
        title('CPU time in function evaluations')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); xlim([0.5,n+0.5]);
        t1_feval(i,j) = t_feval_1ij;
        t2_feval(i,j) = t_feval_2ij;


        subplot(nv,nh,5)
        plot(i,R3ij.stats.t_wall_total-t_feval_1ij,'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
        hold on
        plot(i,R2ij.stats.t_wall_total-t_feval_2ij,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('trial')
        ylabel('time (s)')
        title('CPU time in IPOPT')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); xlim([0.5,n+0.5]);


        subplot(nv,nh,6)
        plot(i,t_feval_1ij/R3ij.stats.t_wall_total*100,'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
        hold on
        plot(i,t_feval_2ij/R2ij.stats.t_wall_total*100,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('trial')
        ylabel('relative time (%)')
        title({'Fraction of CPU time in','function evaluations'})
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); xlim([0.5,n+0.5]);

        subplot(nv,nh,7)
        plot(i,R3ij.stats.iter_count,'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
        hold on
        plot(i,R2ij.stats.iter_count,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('trial')
        ylabel('# iterations')
        title('Number of iterations')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); xlim([0.5,n+0.5]);
        t1_iter(i,j) = R3ij.stats.iter_count;
        t2_iter(i,j) = R2ij.stats.iter_count;

        t1_rel(i,j) = (R3ij.stats.t_wall_total)/R3ij.stats.iter_count;
        t2_rel(i,j) = (R2ij.stats.t_wall_total)/R2ij.stats.iter_count;

        subplot(nv,nh,8)
        plot(i,R3ij.stats.n_call_nlp_f,'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
        hold on
        plot(i,R2ij.stats.n_call_nlp_f,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('trial')
        ylabel('# obj fcn eval')
        title('# Objective function evaluations')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); xlim([0.5,n+0.5]);
        t1_f(i,j) = R3ij.stats.n_call_nlp_f;
        t2_f(i,j) = R2ij.stats.n_call_nlp_f;

        subplot(nv,nh,9)
        plot(i,R3ij.stats.n_call_nlp_grad_f,'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
        hold on
        plot(i,R2ij.stats.n_call_nlp_grad_f,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('trial')
        ylabel('# obj grad eval')
        title('# Objective gradient evaluations')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); xlim([0.5,n+0.5]);
        t1_gf(i,j) = R3ij.stats.n_call_nlp_grad_f;
        t2_gf(i,j) = R2ij.stats.n_call_nlp_grad_f;

        subplot(nv,nh,10)
        plot(i-0.1,R3ij.stats.iterations.obj(end),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
        hold on
        plot(i+0.1,R2ij.stats.iterations.obj(end),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('trial')
        ylabel('objective')
        title('Objective value')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); xlim([0.5,n+0.5]);
        t1_obj(i,j) = R3ij.stats.iterations.obj(end);
        t2_obj(i,j) = R2ij.stats.iterations.obj(end);

        subplot(nv,nh,11)
        rms_obj = rms(R3ij.stats.iterations.obj(end) - R2ij.stats.iterations.obj(end));
        semilogy(i,rms_obj,'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
        hold on
        xlabel('trial')
        ylabel('rmse objective')
        title('Objective value RMS')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)/10,yl(2)*10]); xlim([0.5,n+0.5]);

    end
    subplot(nv,nh,2)
    plot(i+0.2,t1_tot_i,'o','Color',cs3(end,:),'MarkerFaceColor',cs3(end,:))
    hold on
    plot(i+0.2,t2_tot_i,'d','Color',cs2(end,:),'MarkerFaceColor',cs2(end,:))
    xlabel('trial')
    ylabel('time (s)')
    title('Summed CPU time')
    axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); xlim([0.5,n+0.5]);

end


subplot(nv,nh,12)
for j=1:length(noise_factors)
    hold on
    p1=plot(1,noise_factors(j),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:));
    p2=plot(2,noise_factors(j),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:));
    p3=plot(3,noise_factors(j),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:));
end
ylabel('noise factor')
legend([p1,p2,p3],{'further adapted code','adapted code','difference'},'Location','southoutside')
title('Legend')
xlim([0.5,3.5])
ylim([0.5,7.5])

%%

figure('Position',[0,40,scs(3),scs(4)-120]);

for j=1:length(noise_factors)
    subplot(2,5,1)
    plot(noise_factors(j),mean(t1_tot(:,j)),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
    hold on
    plot(noise_factors(j),mean(t2_tot(:,j)),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
    xlabel('noise factor')
    ylabel('time (s)')
    title({'Mean CPU time','total'})
    xlim([0.5,7.5])

    subplot(2,5,2)
    plot(noise_factors(j),mean(t1_feval(:,j)),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
    hold on
    plot(noise_factors(j),mean(t2_feval(:,j)),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
    xlabel('noise factor')
    ylabel('time (s)')
    title({'Mean CPU time','function evaluations'})
    xlim([0.5,7.5])

    subplot(2,5,3)
    plot(noise_factors(j),mean(t1_tot(:,j)-t1_feval(:,j)),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
    hold on
    plot(noise_factors(j),mean(t2_tot(:,j)-t2_feval(:,j)),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
    xlabel('noise factor')
    ylabel('time (s)')
    title({'Mean CPU time','IPOPT'})
    xlim([0.5,7.5])

    subplot(2,5,4)
    plot(noise_factors(j),mean(t1_feval(:,j)./t1_tot(:,j)*100),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
    hold on
    plot(noise_factors(j),mean(t2_feval(:,j)./t2_tot(:,j)*100),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
    xlabel('noise factor')
    ylabel('relative time (%)')
    title({'Mean fraction of CPU time','in function evaluations'})
    xlim([0.5,7.5])

    subplot(2,5,5)
    plot(noise_factors(j),mean(t1_iter(:,j)),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
    hold on
    plot(noise_factors(j),mean(t2_iter(:,j)),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
    xlabel('noise factor')
    ylabel('# iterations')
    title('Number of iterations')
    xlim([0.5,7.5])

    subplot(2,5,6)
    plot(noise_factors(j),mean(t1_f(:,j)),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
    hold on
    plot(noise_factors(j),mean(t2_f(:,j)),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
    xlabel('noise factor')
    ylabel('# obj fcn eval')
    title({'Number of evaluations','objective function'})
    xlim([0.5,7.5])

    subplot(2,5,7)
    plot(noise_factors(j),mean(t1_gf(:,j)),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
    hold on
    plot(noise_factors(j),mean(t2_gf(:,j)),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
    xlabel('noise factor')
    ylabel('# obj grad eval')
    title({'Number of evaluations','objective gradient'})
    xlim([0.5,7.5])

    subplot(2,5,8)
    plot(noise_factors(j),mean(t1_obj(:,j)),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
    hold on
    plot(noise_factors(j),mean(t2_obj(:,j)),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
    xlabel('noise factor')
    ylabel('objective')
    title('Objective value')
    xlim([0.5,7.5])

    subplot(2,5,9)
    plot(noise_factors(j),mean(t1_rel(:,j)),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:))
    hold on
    plot(noise_factors(j),mean(t2_rel(:,j)),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
    xlabel('noise factor')
    ylabel('time / iter (s)')
    title({'Mean CPU time','per iteration'})
    xlim([0.5,7.5])

    subplot(2,5,10)
    p1=plot(1,noise_factors(j),'o','Color',cs3(j,:),'MarkerFaceColor',cs3(j,:));
    hold on
    p2=plot(2,noise_factors(j),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:));
    title('Legend')
    ylabel('noise factor')
    xlim([0.5,2.5])
    ylim([0.5,7.5])

end

subplot(2,5,10)
legend([p1,p2],{'further adapted code','adapted code'},'Location','southoutside')


%% plot results

R31 = Results3{1};
R21 = Results2{1};
dt = 0.167; T = 10.02; t = 0:dt:T;

figure('Position',[0,40,scs(3)/2,scs(4)-120]);

for j=1:length(noise_factors)
    
    R31j = R31.(['R' num2str(noise_factors(j))]);
    R21j = R21.(['R' num2str(noise_factors(j))]);

    subplot(3,2,1)
    plot(t,R31j.x_mean(1,:),'Color',cs3(j,:));
    title('position');
    hold on;
    plot(t,R21j.x_mean(1,:),'--','Color',cs2(j,:));
    xlim([t(1),t(end)])

    subplot(3,2,2)
    plot(t,180/pi*sqrt(R31j.P(1,:)),'Color',cs3(j,:));
    title('position SD');
    hold on;
    plot(t,180/pi*sqrt(R21j.P(1,:)),'--','Color',cs2(j,:));
    set(gca, 'YScale', 'log')
    xlim([t(1),t(end)])

    subplot(3,2,3)
    plot(t,R31j.x_mean(2,:),'Color',cs3(j,:));
    title('velocity');
    hold on;
    plot(t,R21j.x_mean(2,:),'--','Color',cs2(j,:));
    xlim([t(1),t(end)])
    
    subplot(3,2,4)
    plot(t,180/pi*sqrt(R31j.P(4,:)),'Color',cs3(j,:));
    title('velocity SD');
    hold on;
    plot(t,180/pi*sqrt(R21j.P(4,:)),'--','Color',cs2(j,:));
    set(gca, 'YScale', 'log')
    xlim([t(1),t(end)])
    
    subplot(3,2,5)
    plot(t,R31j.F(1,:),'Color',cs3(j,:));
    title('torque');
    hold on;
    plot(t,R21j.F(1,:),'--','Color',cs2(j,:));
    xlim([t(1),t(end)]) 

    subplot(3,2,6)
    plot([1,2],[1,1]*noise_factors(j),'Color',cs3(j,:));
    hold on
    plot([2,3],[1,1]*noise_factors(j),'--','Color',cs2(j,:));
    title('legend');
    ylabel('noise factor')
    xlabel('further adapted    adapted   ')
    ylim([0.5,7.5])
    xlim([0.8,3.2])
    set(gca,'XTickLabel','')

end

%% plot results difference

figure('Position',[scs(3)/2,40,scs(3)/2,scs(4)-120]);

R31 = Results3{1};
R21 = Results2{1};

for j=1:length(noise_factors)
    R31j = R31.(['R' num2str(noise_factors(j))]);
    R21j = R21.(['R' num2str(noise_factors(j))]);

    subplot(3,2,1)
    semilogy(j,rms(R2ij.x_mean(1,:) - R3ij.x_mean(1,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
    hold on
    xlabel('noise factor')
    ylabel('rmse x_m')
    title('RMSE mean position')
    axis tight
    
    subplot(3,2,3)
    semilogy(j,rms(R2ij.x_mean(2,:) - R3ij.x_mean(2,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
    hold on
    xlabel('noise factor')
    ylabel('rmse xdot_m')
    title('RMSE mean velocity')
    
    subplot(3,2,2)
    semilogy(j,rms(R2ij.P(1,:) - R3ij.P(1,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
    hold on
    xlabel('noise factor')
    ylabel('rmse P(1,1)')
    title('RMSE position variance')
    
    subplot(3,2,4)
    semilogy(j,rms(R2ij.P(4,:) - R3ij.P(4,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
    hold on
    xlabel('noise factor')
    ylabel('rmse P(2,2)')
    title('RMSE velocity variance')
    
    subplot(3,2,6)
    semilogy(j,rms(R2ij.P(2,:) - R3ij.P(2,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
    hold on
    xlabel('noise factor')
    ylabel('rmse P(1,2)')
    title('RMSE covariance')
    
    subplot(3,2,5)
    semilogy(j,rms(R2ij.F - R3ij.F),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
    hold on
    xlabel('noise factor')
    ylabel('rmse F')
    title('RMSE control')
end





