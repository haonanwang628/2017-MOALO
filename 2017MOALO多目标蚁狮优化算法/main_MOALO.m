%_________________________________________________________________________%
%  Multi-Objective Ant Lion Optimizer (MALO) source codes demo            %
%                           version 1.0                                   %
%                                                                         %
%  Developed in MATLAB R2011b(7.13)                                       %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
% Paper: Mirjalili, Seyedali, Pradeep Jangir, and Shahrzad Saremi.        %
% Multi-objective ant lion optimizer: a multi-objective optimization      %
%  algorithm for solving engineering problems." Applied Intelligence      %
% (2016): 1-17, DOI: http://dx.doi.org/10.1007/s10489-016-0825-8          %
%_________________________________________________________________________%

% clc;
% clear;
% close all;
% 
% % Change these details with respect to your problem%%%%%%%%%%%%%%
% ObjectiveFunction=@ZDT1;
% dim=5;
% lb=0;
% ub=1;
% obj_no=2;
% 
% if size(ub,2)==1
%     ub=ones(1,dim)*ub;
%     lb=ones(1,dim)*lb;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Initial parameters of the MODA algorithm
% max_iter=100;
% N=100;
% ArchiveMaxSize=100;

clear all;
close all;
clc;

%%%% ����ʵ�������Χ������ʵ�������ʵ�麯����Χ
Num_Test=5;   %%%% ÿ��������������Num_Test��?
Num_Experiment=30;   %%%% �����Ǵ�F1-FNum_Functions
AlgorithmName='MOALO'; %%% ���ƺ�����?


% gmax = 100;    %����������
% % FEE = 5000; %%%���Ŀ�꺯�����۴���
% n = 50;       %��Ⱥ��ģ
max_iter=200;  % Maximum Number of Iterations
N=200;    % Population Size (Number of Sub-Problems)
ArchiveMaxSize=100;  %%%Archive Size(number of rep)

m = 5;   %Ŀ��ά��

ALLFunction_AllTest=[];

% for ff=[1:21];
% for ff=[10:21];
% for ff=[1:4,6:9];
for ff=[7];
    clearvars -except Num_Test Num_Experiment AlgorithmName ALLFunction_AllTest ff max_iter N ArchiveMaxSize m AllTest_Results  problem_name 
    %%%%% �����ļ���?
    string_0ALL=['000\',AlgorithmName,'_5άĿ��800�ε���100��Ⱥʵ��20210923\'];
    dirname00=[string_0ALL,'\F',num2str(ff),'\'];
   display(['**********  ',AlgorithmName,'�㷨�Ż�F',num2str(ff),'�� ', 'M',num2str(m), ' άʵ��   **********']);
   for testi=1:Num_Test   %%%% ����ÿ��ʵ����Դ���
       dirname0=[dirname00,'test',num2str(testi),'_F',num2str(ff)];
       system(['mkdir ' dirname0]) %�������ļ���
       dirname1=[dirname0,'\F',num2str(ff),'_fig'];
       system(['mkdir ' dirname1]) %�����ļ���  �ȴ�����ʵ��ͼ��
       dirname2=[dirname0,'\F',num2str(ff),'_data'];
       system(['mkdir ' dirname2]) %�����ļ���  �ȴ�����ʵ��ͼ��
       for kk=1:30 %%%% ����ʵ�������ѭ��
           display(['**********  ',AlgorithmName,'�㷨�Ż�F',num2str(ff),'��  ��  ', num2str(kk), ' ��ʵ��   **********']);
            rand('state',sum(100*clock));
            problem_name=['F',num2str(ff)];
   
           [ ub,lb,dim ] = generate_boundary1( problem_name,m );%Upper and Lower Bound of Decision Variables  %%%���ɾ��߿ռ��б����Ͻ硢�½��ά��
tic;  % CPU time measure
Archive_X=zeros(100,dim);
Archive_F=ones(100,m)*inf;

Archive_member_no=0;

r=(ub-lb)/2;
V_max=(ub(1)-lb(1))/10;

Elite_fitness=inf*ones(1,m);
Elite_position=zeros(dim,1);

Ant_Position=initialization(N,dim,ub,lb);
fitness=zeros(N,2);

V=initialization(N,dim,ub,lb);
iter=0;

position_history=zeros(N,max_iter,dim);

for iter=1:max_iter
    
    for i=1:N %Calculate all the objective values first
%         Particles_F(i,:)=ObjectiveFunction(Ant_Position(:,i)');
        Particles_F(i,:)=test_function(Ant_Position(:,i)',dim,m,problem_name);
        if dominates(Particles_F(i,:),Elite_fitness)
            Elite_fitness=Particles_F(i,:);
            Elite_position=Ant_Position(:,i);
        end
    end
    
    [Archive_X, Archive_F, Archive_member_no]=UpdateArchive(Archive_X, Archive_F, Ant_Position, Particles_F, Archive_member_no);
    
    if Archive_member_no>ArchiveMaxSize
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, m);
        [Archive_X, Archive_F, Archive_mem_ranks, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
    else
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, m);
    end
    
    Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, m);
    
    % Chose the archive member in the least population area as arrtactor
    % to improve coverage
    index=RouletteWheelSelection(1./Archive_mem_ranks);
    if index==-1
        index=1;
    end
    Elite_fitness=Archive_F(index,:);
    Elite_position=Archive_X(index,:)';
    
    Random_antlion_fitness=Archive_F(1,:);
    Random_antlion_position=Archive_X(1,:)';
    
    for i=1:N
        
        index=0;
        neighbours_no=0;
        
        RA=Random_walk_around_antlion(dim,max_iter,lb,ub, Random_antlion_position',iter);
        
        [RE]=Random_walk_around_antlion(dim,max_iter,lb,ub, Elite_position',iter);
        
        Ant_Position(:,i)=(RE(iter,:)'+RA(iter,:)')/2;
        
        
        
        Flag4ub=Ant_Position(:,i)>ub';
        Flag4lb=Ant_Position(:,i)<lb';
        Ant_Position(:,i)=(Ant_Position(:,i).*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
        
    end
%     display(['At the iteration ', num2str(iter), ' there are ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);
    HisPF{iter} = Archive_F;
end

% figure
% 
% Draw_ZDT1();
% 
% hold on
% 
% plot(Archive_F(:,1),Archive_F(:,2),'ko','MarkerSize',8,'markerfacecolor','k');
% 
% legend('True PF','Obtained PF');
% title('MALO');
% 
% set(gcf, 'pos', [403   466   230   200])

            time=toc;
            PF = Archive_F;
            cg_curve=HisPF; %%% ��ʷĿ�꺯��ֵ
            Time(kk)=time;
            cc=strcat(dirname2,'\',AlgorithmName,'�Ż�����_',num2str(kk),'.mat');
            result.time=time;
           
            true_PF=TPF(m,ArchiveMaxSize, problem_name);

         %%% using the matlab codes for calculating metric values
         

            hv = HV(PF,true_PF);   %�����?
            gd = GD(PF, true_PF);                       %�������릻

            sp = Spacing(PF, true_PF);                  %�ռ�ֲ� 
            igd = IGD(PF, true_PF);            %�����������릻

            hvd(kk)=hv;
            gdd(kk)=gd;
            ssp(kk)=sp;
            igdd(kk)=igd;

             save(cc)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       end
       mean_IGD = mean(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���IGDƽ��ֵ : ', num2str(mean_IGD)]);
       std_IGD=std(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���IGD��׼�� : ', num2str(std_IGD)]);
       max_IGD=max(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���IGD���ֵ : ', num2str(max_IGD)]);
       min_IGD=min(igdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���IGD��Сֵ : ', num2str(min_IGD)]);
       display('******************************** ');
       
       
       mean_GD = mean(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���GDƽ��ֵ : ', num2str(mean_GD)]);
       std_GD=std(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���GD��׼�� : ', num2str(std_GD)]);
       max_GD=max(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���GD���ֵ : ', num2str(max_GD)]);
       min_GD=min(gdd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���GD��Сֵ : ', num2str(min_GD)]);
       display('******************************** ');
      
       
       mean_HV = mean(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���HVƽ��ֵ : ', num2str(mean_HV)]);
       std_HV=std(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���HV��׼�� : ', num2str(std_HV)]);
       max_HV=max(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���HV���ֵ : ', num2str(max_HV)]);
       min_HV=min(hvd);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���HV��Сֵ : ', num2str(min_HV)]);
       display('******************************** ');
       
       mean_SP = mean(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���SPƽ��ֵ : ', num2str(mean_SP)]);
       std_SP=std(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���SP��׼�� : ', num2str(std_SP)]);
       max_SP=max(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���SP���ֵ : ', num2str(max_SP)]);
       min_SP=min(ssp);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ���SP��Сֵ : ', num2str(min_SP)]);
       display('******************************** ');
       
       mean_time=mean(Time);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ�������ʱ��ƽ��ֵ : ', num2str(mean_time)]);
       std_time=std(Time);
       display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ�������ʱ���׼�� : ', num2str(std_time)]);
       display('******************************** ');
%        mean_X=mean(Best_X);
%         display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ������Ž�ƽ��ֵ : ', num2str(mean_X)]);
%         std_X=std(Best_X);
%         display([AlgorithmName,'_Functions_F',num2str(ff),'����',num2str(kk),'��ʵ������Ž��׼�� : ', num2str(std_X)]);
%         %%%%%%%%%%%%%%%%%%
        cd=strcat(dirname0,'\Result���ܽ��.mat');
        Result.IGDmean=mean_IGD;
        Result.IGDstd=std_IGD;
        Result.IGDmax=max_IGD;
        Result.IGDmin=min_IGD;
      
        
        Result.GDmean=mean_GD;
        Result.GDstd=std_GD;
        Result.GDmax=max_GD;
        Result.GDmin=min_GD;
        
        Result.HVmean=mean_HV;
        Result.HVstd=std_HV;
        Result.HVmax=max_HV;
        Result.HVmin=min_HV;
        
        Result.SPmean=mean_SP;
        Result.SPstd=std_SP;
        Result.SPmax=max_SP;
        Result.SPmin=min_SP;
        
        Result.Tmean=mean_time;
        Result.Tstd=std_time;
        %         Result.Xmean=mean_X;
        %         Result.Xstd=std_X;
        %         Result.Best_Y=Best_Y;
        %         Result.Best_X=Best_X;
        Result.Time=Time;
        %         Result.ResultVector=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_GD,std_GD,max_GD,min_GD,mean_time,std_time];
        Result.ResultVector=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_GD,std_GD,max_GD,min_GD,mean_HV,std_HV,max_HV,min_HV,mean_SP,std_SP,max_SP,min_SP,mean_time,std_time];
        %         Result.Best_History_Y=History_Y;
        save(cd,'Result')
        
        
        %         AllTest_Results(testi,:)=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_time,std_time];
        % AllTest_Results(testi,:)=[mean_IGD,std_IGD,max_IGD,min_IGD,mean_GD,std_GD,max_GD,min_GD,mean_time,std_time];
        AllTest_Results(testi,:)=[mean_IGD,std_IGD,mean_GD,std_GD,mean_HV,std_HV,mean_SP,std_SP,mean_time,std_time];
    end
    cd=strcat(dirname00,'Result_AllTest.mat');
    save(cd,'AllTest_Results')
    ALLFunction_AllTest=[ALLFunction_AllTest;AllTest_Results];
end


cd=strcat(string_0ALL,'ALLFunction_AllTest.mat');
save(cd,'ALLFunction_AllTest')