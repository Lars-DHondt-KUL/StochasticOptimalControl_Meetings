clear
% close all
clc

addpath([pwd '\Integrator']);



noise_factors = [2];
dt = 0.167;
dt_sf = [0.5,1,1.5,2];
n = length(dt_sf);
res_name = ['compare_noise_factor_' num2str(noise_factors(1)) '_trials_' num2str(n), '_dt'];


%%
% time_1 = nan(n,1);
% time_2 = nan(n,1);
% 
% Results1 = cell(n,1);
% Results2 = cell(n,1);

% log_folder = ['logfiles_' res_name];
% mkdir(log_folder);

%%
% for i=1:n
% 
%     dti = dt/dt_sf(i);
%     
%     diary([pwd '/' log_folder '/' 'example2_2_i' num2str(i) '.txt'])
%     t0 = tic;
%     [R1] = example2_2(0,noise_factors,dti);
%     t1 = toc(t0);
%     Results1{i} = R1;
%     time_1(i) = t1;
%     diary off
% 
% 
%     diary([pwd '/' log_folder '/' 'example2_2b_i' num2str(i) '.txt'])
%     t0 = tic;
%     [R2] = example2_2b(0,noise_factors,dti);
%     t2 = toc(t0);
%     Results2{i} = R2;
%     time_2(i) = t2;
%     diary off
% 
% end

%%

% save([res_name '.mat'],'time_1','time_2','Results1','Results2');

load([res_name '.mat'])

%%

%% plot comparison

cs1 = linspecer(length(noise_factors)+1,'blue');
cs1 = cs1(2:end,:);

cs2 = linspecer(length(noise_factors)+1,'green');
cs2 = cs2(2:end,:);

csd = linspecer(length(noise_factors)+1,'red');
csd = csd(2:end,:);

nv = 2;
nh = 7;

scs = get(0,'ScreenSize');

figure('Position',[0,40,scs(3),scs(4)-120]);



for i=1:n
    R1i = Results1{i};
    R2i = Results2{i};

    t1_tot_i =0;
    t2_tot_i =0;

    for j=1:length(noise_factors)
        R1ij = R1i.(['R' num2str(noise_factors(j))]);
        R2ij = R2i.(['R' num2str(noise_factors(j))]);

        xn = length(R2ij.F)-1;

        if j==1
            subplot(nv,nh,1)
            plot(xn,time_1(i),'o','Color',cs1(end,:),'MarkerFaceColor',cs1(end,:))
            hold on
            plot(xn,time_2(i),'d','Color',cs2(end,:),'MarkerFaceColor',cs2(end,:))
            xlabel('N time intervals')
            ylabel('time (s)')
            title('Total trial duration')
            axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
            xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
            set(gca,'XTick',linspace(xl(1),xl(2),n))
        end

        subplot(nv,nh,3)
        plot(xn,R1ij.stats.t_wall_total,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,R2ij.stats.t_wall_total,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('time (s)')
        title('Total CPU time')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))

        t1_tot_i = t1_tot_i + R1ij.stats.t_wall_total;
        t2_tot_i = t2_tot_i + R2ij.stats.t_wall_total;


        subplot(nv,nh,4)
        t_feval_1ij = R1ij.stats.t_wall_nlp_f + R1ij.stats.t_wall_nlp_g + ...
            R1ij.stats.t_wall_nlp_grad + R1ij.stats.t_wall_nlp_grad_f + ...
            R1ij.stats.t_wall_nlp_hess_l + R1ij.stats.t_wall_nlp_jac_g;
        t_feval_2ij = R2ij.stats.t_wall_nlp_f + R2ij.stats.t_wall_nlp_g + ...
            R2ij.stats.t_wall_nlp_grad + R2ij.stats.t_wall_nlp_grad_f + ...
            R2ij.stats.t_wall_nlp_hess_l + R2ij.stats.t_wall_nlp_jac_g;
        plot(xn,t_feval_1ij,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,t_feval_2ij,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('time (s)')
        title('CPU time in function evaluations')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))


        subplot(nv,nh,5)
        plot(xn,R1ij.stats.t_wall_total-t_feval_1ij,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,R2ij.stats.t_wall_total-t_feval_2ij,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('time (s)')
        title('CPU time in IPOPT')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))


        subplot(nv,nh,6)
        plot(xn,t_feval_1ij/R1ij.stats.t_wall_total*100,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,t_feval_2ij/R2ij.stats.t_wall_total*100,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('relative time (%)')
        title({'Fraction of CPU time in','function evaluations'})
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))

        subplot(nv,nh,7)
        plot(xn,R1ij.stats.iter_count,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,R2ij.stats.iter_count,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('# iterations')
        title('Number of iterations')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))

        subplot(nv,nh,8)
        plot(xn,R1ij.stats.n_call_nlp_f,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,R2ij.stats.n_call_nlp_f,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('# obj fcn eval')
        title('# Objective function evaluations')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))

        subplot(nv,nh,9)
        plot(xn,R1ij.stats.n_call_nlp_grad_f,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,R2ij.stats.n_call_nlp_grad_f,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('# obj grad eval')
        title('# Objective gradient evaluations')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))

        subplot(nv,nh,10)
        plot(xn,R1ij.stats.iterations.obj(end),'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,R2ij.stats.iterations.obj(end),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('objective')
        title('Objective value')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))

        subplot(nv,nh,11)
        rms_obj = rms(R1ij.stats.iterations.obj(end) - R2ij.stats.iterations.obj(end));
        semilogy(xn,rms_obj,'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))%
        hold on
        xlabel('N time intervals')
        ylabel('rmse objective')
        title('Objective value RMS')
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)/10,yl(2)*10]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))

        subplot(nv,nh,12)
        t_feval_1ij = (R1ij.stats.t_wall_nlp_f + R1ij.stats.t_wall_nlp_g)/R1ij.stats.n_call_nlp_f;
        t_feval_2ij = (R2ij.stats.t_wall_nlp_f + R2ij.stats.t_wall_nlp_g)/R2ij.stats.n_call_nlp_f;
            
        plot(xn,t_feval_1ij*1e3,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,t_feval_2ij*1e3,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('time / nlp eval (ms)')
        title({'Average CPU time', 'per nlp evaluation'})
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))

        subplot(nv,nh,13)
        t_feval_1ij = (R1ij.stats.t_wall_total)/R1ij.stats.iter_count;
        t_feval_2ij = (R2ij.stats.t_wall_total)/R2ij.stats.iter_count;
            
        plot(xn,t_feval_1ij,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
        hold on
        plot(xn,t_feval_2ij,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
        xlabel('N time intervals')
        ylabel('time / iter (s)')
        title({'Average CPU time','per iteration'})
        axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
        xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
        set(gca,'XTick',linspace(xl(1),xl(2),n))


    end
    subplot(nv,nh,2)
    plot(xn,t1_tot_i,'o','Color',cs1(end,:),'MarkerFaceColor',cs1(end,:))
    hold on
    plot(xn,t2_tot_i,'d','Color',cs2(end,:),'MarkerFaceColor',cs2(end,:))
    xlabel('N time intervals')
    ylabel('time (s)')
    title('Summed CPU time')
    axis tight; yl = get(gca, 'ylim'); ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)]); 
    xl = get(gca, 'xlim'); xlim([xl(1)-0.1*norm(xl),xl(2)+0.1*norm(xl)]);
    set(gca,'XTick',linspace(xl(1),xl(2),n))

end


subplot(nv,nh,nv*nh)
for j=1:length(noise_factors)
    hold on
    p1=plot(1,noise_factors(j),'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:));
    p2=plot(2,noise_factors(j),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:));
    p3=plot(3,noise_factors(j),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:));
end
ylabel('noise factor')
legend([p1,p2,p3],{'original code','adapted code','difference'},'Location','southoutside')
title('Legend')
xlim([0.5,3.5])
ylim([noise_factors(1)-0.5,noise_factors(1)+0.5])


%%
% 
% cs1 = linspecer(length(noise_factors)+1,'blue');
% cs1 = cs1(2:end,:);
% 
% cs2 = linspecer(length(noise_factors)+1,'green');
% cs2 = cs2(2:end,:);
% 
% csd = linspecer(length(noise_factors)+1,'red');
% csd = csd(2:end,:);
% 
% nv = 3;
% nh = 5;
% 
% scs = get(0,'ScreenSize');
% 
% figure('Position',[0,40,scs(3),scs(4)-140]);
% 
% subplot(nv,nh,1)
% pl1=plot(1:n,time_1,'o','Color',cs1(end,:),'MarkerFaceColor',cs1(end,:),'DisplayName','elapsed time');
% hold on
% pl2=plot(1:n,time_2,'d','Color',cs2(end,:),'MarkerFaceColor',cs2(end,:),'DisplayName','elapsed time');
% xlabel('N time intervals')
% ylabel('time (s)')
% title('Total duration')
% xlim([0.5,n+0.5])
% 
% for i=1:n
%     R1i = Results1{i};
%     R2i = Results2{i};
% 
%     t1_tot_i =0;
%     t2_tot_i =0;
% 
%     for j=1:length(noise_factors)
%         R1ij = R1i.(['R' num2str(noise_factors(j))]);
%         R2ij = R2i.(['R' num2str(noise_factors(j))]);
% 
%         subplot(nv,nh,2)
%         plot(i,R1ij.stats.t_wall_total,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
%         hold on
%         plot(i,R2ij.stats.t_wall_total,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
%         xlabel('N time intervals')
%         ylabel('time (s)')
%         title('Total CPU time')
%         xlim([0.5,n+0.5])
% 
%         t1_tot_i = t1_tot_i + R1ij.stats.t_wall_total;
%         t2_tot_i = t2_tot_i + R2ij.stats.t_wall_total;
% 
%         subplot(nv,nh,3)
%         t_feval_1ij = R1ij.stats.t_wall_nlp_f + R1ij.stats.t_wall_nlp_g + ...
%             R1ij.stats.t_wall_nlp_grad + R1ij.stats.t_wall_nlp_grad_f + ...
%             R1ij.stats.t_wall_nlp_hess_l + R1ij.stats.t_wall_nlp_jac_g;
%         t_feval_2ij = R2ij.stats.t_wall_nlp_f + R2ij.stats.t_wall_nlp_g + ...
%             R2ij.stats.t_wall_nlp_grad + R2ij.stats.t_wall_nlp_grad_f + ...
%             R2ij.stats.t_wall_nlp_hess_l + R2ij.stats.t_wall_nlp_jac_g;
%         plot(i,t_feval_1ij,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
%         hold on
%         plot(i,t_feval_2ij,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
%         xlabel('N time intervals')
%         ylabel('time (s)')
%         title('Total CPU time in function evaluations')
%         xlim([0.5,n+0.5])
% 
%         subplot(nv,nh,4)
%         plot(i,R1ij.stats.t_wall_total-t_feval_1ij,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
%         hold on
%         plot(i,R2ij.stats.t_wall_total-t_feval_2ij,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
%         xlabel('N time intervals')
%         ylabel('time (s)')
%         title('Total CPU time in IPOPT')
%         xlim([0.5,n+0.5])
% 
%         subplot(nv,nh,5)
%         plot(i,t_feval_1ij/R1ij.stats.t_wall_total*100,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
%         hold on
%         plot(i,t_feval_2ij/R2ij.stats.t_wall_total*100,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
%         xlabel('N time intervals')
%         ylabel('relative time (%)')
%         title('Fraction of CPU time in function evaluations')
%         xlim([0.5,n+0.5])
% 
%         subplot(nv,nh,6)
%         plot(i,R1ij.stats.iter_count,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
%         hold on
%         plot(i,R2ij.stats.iter_count,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
%         xlabel('N time intervals')
%         ylabel('# iterations')
%         title('Number of iterations')
%         xlim([0.5,n+0.5])
% 
%         subplot(nv,nh,7)
%         plot(i,R1ij.stats.n_call_nlp_f,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
%         hold on
%         plot(i,R2ij.stats.n_call_nlp_f,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
%         xlabel('N time intervals')
%         ylabel('# obj fcn eval')
%         title('# Objective function evaluations')
%         xlim([0.5,n+0.5])
% 
%         subplot(nv,nh,8)
%         plot(i,R1ij.stats.n_call_nlp_grad_f,'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
%         hold on
%         plot(i,R2ij.stats.n_call_nlp_grad_f,'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
%         xlabel('N time intervals')
%         ylabel('# obj grad eval')
%         title('# Objective gradient evaluations')
%         xlim([0.5,n+0.5])
% 
%         subplot(nv,nh,9)
%         plot(i-0.1,R1ij.stats.iterations.obj(end),'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:))
%         hold on
%         plot(i+0.1,R2ij.stats.iterations.obj(end),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:))
%         xlabel('N time intervals')
%         ylabel('objective')
%         title('Objective value')
%         xlim([0.5,n+0.5])
% 
%         subplot(nv,nh,10)
%         rms_obj = rms(R1ij.stats.iterations.obj(end) - R2ij.stats.iterations.obj(end));
%         semilogy(i,rms_obj,'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
%         hold on
%         xlabel('N time intervals')
%         ylabel('objective')
%         title('Objective value')
%         xlim([0.5,n+0.5])
% 
%         if i==1
%             subplot(nv,6,(nv-1)*6+1)
%             semilogy(j,rms(R2ij.x_mean(1,:) - R1ij.x_mean(1,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
%             hold on
%             xlabel('noise factor')
%             ylabel('rmse x_m')
%             title('RMSE mean position')
%     
%             subplot(nv,6,(nv-1)*6+2)
%             semilogy(j,rms(R2ij.x_mean(2,:) - R1ij.x_mean(2,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
%             hold on
%             xlabel('noise factor')
%             ylabel('rmse xdot_m')
%             title('RMSE mean velocity')
%     
%             subplot(nv,6,(nv-1)*6+3)
%             semilogy(j,rms(R2ij.P(1,:) - R1ij.P(1,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
%             hold on
%             xlabel('noise factor')
%             ylabel('rmse P(1,1)')
%             title('RMSE position variance')
%     
%             subplot(nv,6,(nv-1)*6+4)
%             semilogy(j,rms(R2ij.P(4,:) - R1ij.P(4,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
%             hold on
%             xlabel('noise factor')
%             ylabel('rmse P(2,2)')
%             title('RMSE velocity variance')
% 
%             subplot(nv,6,(nv-1)*6+5)
%             semilogy(j,rms(R2ij.P(2,:) - R1ij.P(2,:)),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
%             hold on
%             xlabel('noise factor')
%             ylabel('rmse P(1,2)')
%             title('RMSE covariance')
% 
%             subplot(nv,6,(nv-1)*6+6)
%             semilogy(j,rms(R2ij.F - R1ij.F),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:))
%             hold on
%             xlabel('noise factor')
%             ylabel('rmse F')
%             title('RMSE control')
% 
%         end
% 
%     end
%     subplot(nv,nh,1)
%     pl3=plot(i+0.2,t1_tot_i,'o','Color',cs1(end,:),'DisplayName','summed CPU time');
%     pl4=plot(i+0.2,t2_tot_i,'d','Color',cs2(end,:),'DisplayName','summed CPU time');
% 
% end
% legend([pl1,pl2,pl3,pl4],'Location','best')

% subplot(nv,nh,10)
% for j=1:length(noise_factors)
%     hold on
%     p1=plot(1,noise_factors(j),'o','Color',cs1(j,:),'MarkerFaceColor',cs1(j,:));
%     p2=plot(2,noise_factors(j),'d','Color',cs2(j,:),'MarkerFaceColor',cs2(j,:));
%     p3=plot(3,noise_factors(j),'v','Color',csd(j,:),'MarkerFaceColor',csd(j,:));
% end
% ylabel('noise factor')
% legend([p1,p2,p3],{'original code','adapted code','difference'},'Location','northeastoutside')
% title('noise factor')


%%
% R11 = Results1{1};
% R21 = Results2{1};
% dt = 0.167; T = 10.02; t = 0:dt:T;
% 
% figure
% for j=1:length(noise_factors)
%     
%     R11j = R11.(['R' num2str(noise_factors(j))]);
%     R21j = R21.(['R' num2str(noise_factors(j))]);
% 
%     subplot(3,2,1)
%     plot(t,R11j.x_mean(1,:),'Color',cs1(j,:));
%     title('position');
%     hold on;
%     plot(t,R21j.x_mean(1,:),'--','Color',cs2(j,:));
%     xlim([t(1),t(end)])
% 
%     subplot(3,2,2)
%     plot(t,180/pi*sqrt(R11j.P(1,:)),'Color',cs1(j,:));
%     title('position SD');
%     hold on;
%     plot(t,180/pi*sqrt(R21j.P(1,:)),'--','Color',cs2(j,:));
%     set(gca, 'YScale', 'log')
%     xlim([t(1),t(end)])
% 
%     subplot(3,2,3)
%     plot(t,R11j.x_mean(2,:),'Color',cs1(j,:));
%     title('velocity');
%     hold on;
%     plot(t,R21j.x_mean(2,:),'--','Color',cs2(j,:));
%     xlim([t(1),t(end)])
%     
%     subplot(3,2,4)
%     plot(t,180/pi*sqrt(R11j.P(4,:)),'Color',cs1(j,:));
%     title('velocity SD');
%     hold on;
%     plot(t,180/pi*sqrt(R21j.P(4,:)),'--','Color',cs2(j,:));
%     set(gca, 'YScale', 'log')
%     xlim([t(1),t(end)])
%     
%     subplot(3,2,5)
%     plot(t,R11j.F(1,:),'Color',cs1(j,:));
%     title('torque');
%     hold on;
%     plot(t,R21j.F(1,:),'--','Color',cs2(j,:));
%     xlim([t(1),t(end)]) 
% end






