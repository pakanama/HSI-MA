% Improved Monkey Algorithm (HSI-MA)

clear
clc

Dim_all = [2 4 6 8 10 15 20];

for Dim = 1:length(Dim_all)
    
    %------Parameter Setting------
    M = 50; % Population size of monkeys
    D = Dim_all(Dim);  % Problem dimension
    MT = 1000; % Allowable iteration number 
    Nc = 100;  % Allowable climb number 
    a = 0.001;  % Step length of climb process
    b = 0.05;  % Eyesight of the monkey
    c = -0.5;  % Lower bound of somersault interval
    d = 0.5;   % Upper bound of somersault interval
    
    %-----Initialization------
    lob(1:D) = -5.12; % Lower population initial boundary
    upb(1:D) = 5.12;  % Upper population initial boundary
    dis=upb-lob;      % Boundary distance
    for iL = 1:1
        % Population initialization
        x=rand(M,D).*repmat(dis,[M,1])+repmat(lob,[M,1]);  
        for kk=1:MT
            %-----Climb Process------
            for i=1:M
                for ite=1:Nc
                    dx(i,:) = (sign(rand(1,D)*2-1))*a;
                    df = (objj(x(i,:)+dx(i,:))-objj(x(i,:)-dx(i,:)))./(2*dx(i,:));  % Pseudo-gradient
                    if ite <= 10
                        Y = x(i,:)-a*sign(df);
                    end
                    if ite > 10
                        Y = x(i,:)-d_tem_dis'.*sign(df);
                    end
                    if sum(lob<Y&Y<upb)==D && objj(Y) < objj(x(i,:))
                        x(i,:) = Y;
                    end
                    % Record location after each climb, and calculate d_ite from the most recent 5 climbs
                    x_ite(ite,i,:) =  x(i,:);
                    if ite>9
                        d_tem_dis = squeeze(max(x_ite(ite-5:ite,i,:))- mean(x_ite(ite-5:ite,i,:))); % Next step length
                    end
                    
                end
            end
            
            %-----Watch Jump Process------
            %-----Repeat until one of the requirements is satisfied: 1. find better location; 2. reach the allowable watch number
            for i=1:M
                if kk>5
                    Max_pos_dis = (max(x_recoder(kk-5:kk-1,i,1:D))-min(x_recoder(kk-5:kk-1,i,1:D)));% Matrix of position change
                    dis_pos_max = max(Max_pos_dis); % Maximum position change
                    dis_obj_max = (objj(x_recoder(kk-5,i,:))-objj(x_recoder(kk-1,i,:)))/objj(x_recoder(kk-5,i,:)); % Objective function value change
                    if dis_pos_max<0.01 && dis_obj_max<0.1 % If one monkey does not change its postion or objective function value much after multiple iterations, send it to Watch
                        % Control the allowable monkeys that do not Watch
                        tcm = 1;
                        tdm = 1;
                        while tcm&&(tdm<5)
                            L_up_tem = lob;
                            L_lo_tem = upb;
                            Y = unifrnd(L_up_tem,L_lo_tem);
                            if objj(Y)<objj(x(i,:))
                                x(i,:) = Y;
                                tcm = 0;
                            end
                            tdm = tdm+1;
                        end
                    end
                end
            end
            
            
            %-----Climb Process------
            for i=1:M
                for ite=1:Nc
                    dx(i,:) = (sign(rand(1,D)*2-1))*a;
                    df = (objj(x(i,:)+dx(i,:))-objj(x(i,:)-dx(i,:)))./(2*dx(i,:));
                    if ite <= 10
                        Y = x(i,:)-a*sign(df);
                    end
                    if ite > 10
                        Y = x(i,:)-d_tem_dis'.*sign(df);
                    end
                    if sum(lob<Y&Y<upb)==D && objj(Y) < objj(x(i,:))
                        x(i,:) = Y;
                    end
                    
                    x_ite(ite,i,:) =  x(i,:);
                    if ite>10
                        d_tem_dis = squeeze(max(x_ite(ite-5:ite,i,:))- mean(x_ite(ite-5:ite,i,:)));
                    end
                end
            end
            
            
            x_recoder(kk,1:M,1:D) = x;
            x_iL_rec_obj(iL,kk) = min(objj(x));
            
            if min(objj(x))<1e-16 ||  kk == MT
                break
            end
            
            %-----Somersault Process------
            tem = 1;
            P = mean(x);   % Population center
            tfm = 1;
            for i=1:M
                if kk>5
                    dis_pos_max = max(max(x_recoder(kk-5:kk-1,i,1:D))-min(x_recoder(kk-5:kk-1,i,1:D))); % Position change
                    dis_obj_max = (objj(x_recoder(kk-5,i,:))-objj(x_recoder(kk-1,i,:)))/objj(x_recoder(kk-5,i,:)); % Objective function value change
                    if dis_pos_max<0.01 &&  dis_obj_max <0.1; 
                        while tem&&(tfm<5)
                            alfa = rand(1,D)*(d-c)+c;  % Control parameter for somersault length control
                            Y = x(i,:)+alfa.*(P-x(i,:));
                            if sum(lob<Y&Y<upb)==D
                                x(i,:) = Y;
                                tem = 0;
                            end
                            tfm=tfm+1;
                        end
                        tem = 1;
                    end
                end
            end
            
            %-----Search Process------
            % Euclidean Distance Matrix
            for im = 1:M
                for ic = 1:M
                    Chrom_dis(im,ic) = norm(x(im,:)-x(ic,:));
                end
            end
            if max(max(Chrom_dis)) < norm(upb-lob)*0.02*M
                x(randi(M),:) = rand(1,D).*repmat(dis,[1,1])+repmat(lob,[1,1]);
            end
            
        end
        
        BestCost_recoder(Dim,iL) = min(min(x_iL_rec_obj(iL,:)'));
    end
    save("Results_"+num2str(D))
    clear x P Y tfm x_recoder x_iL_rec_obj alfa dx dis Chrom_dis dis_obj_max x_ite d_tem_dis
end

