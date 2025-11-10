classdef collocated_segregated
    % collocated_segregated 同位_分离求解器
    properties
        % 求解器状态
        stop_sim = false

    end
    
    methods (Static)
        function [stop_sim, l2_t_over_l2_max_t] = solve(iter, case_obj, fluid, boundary_obj, post)
            % 输入参数:
            %   iter - 迭代次数
            %   case_obj - 案例对象 case是保留字
            %   fluid - 流体属性对象  
            %   boundary_obj - 流体边界对象
            %   post - 后处理对象
            %
            % 输出参数:
            %   stop_sim - 仿真停止标志
            
            arguments
                iter (1,1) {mustBeInteger, mustBePositive}
                case_obj (1,1) structured_mesh
                fluid (1,1) fluid_properties
                boundary_obj (1,1) boundary
                post (1,1) struct
            end
            
            % 存储旧的变量内容
            case_obj.t0 = case_obj.t;
            
            % 从对象获取参数
            dim     = case_obj.dim;
            ncx     = case_obj.ncx;    
            ncy     = case_obj.ncy;
            ncz     = case_obj.ncz;
            ncoef   = case_obj.ncoef;  
            dt      = case_obj.dt;
            
            % 获取场数据
            t       = case_obj.t;
            t0      = case_obj.t0;   
            uf      = case_obj.uf;     
            vf      = case_obj.vf;     
            wf      = case_obj.wf;     
            coef      = case_obj.coef;     
            
            % 求解器参数
            niter_t = case_obj.niter_t;
            relax_t = case_obj.relax_t;
            
            % 方程类型和收敛标准
            ieqn = case_obj.ieqn;
            eqn_conduction = case_obj.eqn_conduction;
            temp_tol = case_obj.temp_tol;
            res_t = case_obj.res_t;
            stop_sim = case_obj.stop_sim;
            
            % 流体属性
            spht     = fluid.spheat;
            con      = fluid.con; 
            dens     = fluid.dens;
            heat_src = fluid.heat_source;
            
            if ieqn == eqn_conduction  % 如果方程类型是温度场类型
                
                % 计算温度方程系数
                coef = collocated_segregated.conduction_coef(case_obj, dim, ncx, ncy, ncz, ncoef, dt, ...
                    spht, con, heat_src, dens, t, uf, vf, wf, coef);
                
                % 处理边界条件
                ct = collocated_segregated.conduction_coef_bound(case_obj, boundary_obj, dim, ncx, ncy, ncz, ...
                    ncoef, dt, con, dens, t, coef);
                
                % 进行雅可比迭代
                initzero = false;
                [t, rel_norm] = solver.scalar_pj(case_obj, post, dim, iter, niter_t, relax_t, ncx, ncy, ncz, ...
                    ncoef, ct, t, initzero, res_t);

                case_obj.t = t;
                case_obj.coef = ct;
                
                % 计算残差范数
                residual_t = t - t0;
                l2_t = norm(residual_t(:), 2);  % 2范数
                if iter > 1
                    l2_max_t = max(l2_t, case_obj.norm2_max_temp);
                    case_obj.norm2_max_temp = l2_max_t;
                else
                    % 第一轮迭代：初始化最大2范数为当前2范数
                    l2_max_t = l2_t;
                    case_obj.norm2_max_temp = l2_t;
                end
                case_obj.norm2_temp = l2_t;
                l2_t_over_l2_max_t = l2_t / l2_max_t;

                % 检查收敛性
                if l2_t_over_l2_max_t < temp_tol && rel_norm > 1 - res_t
                    stop_sim = true;
                    fprintf('\n----------------------------\n');
                    fprintf('Final iter = %d\n', iter);
                    fprintf('it, l2_t/l2_max_t: %d, %.6e\n', iter, l2_t/l2_max_t);
                    fprintf('----------------------------\n');
                end
            end
        end

        function coef_tensor = conduction_coef(case_obj, dim, ncx, ncy, ncz, ncoef, dt, spht, con, ...
                heat_src, dens, t, uf, vf, wf, coef_tensor)
            % 计算温度方程系数
            %
            % 输入参数:
            % case_obj - 案例对象
            % dim - 维度
            % ncx, ncy, ncz - 各方向单元数
            % ncoef - 系数个数
            % dt - 时间步长
            % spht - 比热
            % con - 传导系数
            % heat_src - 热源
            % dens - 密度场
            % t - 温度场
            % uf, vf, wf - 速度场
            % coef_tensor - 系数张量

            arguments
                case_obj (1,1) structured_mesh
                dim (1,1) {mustBeMember(dim, [2, 3])}
                ncx (1,1) {mustBeInteger, mustBePositive}
                ncy (1,1) {mustBeInteger, mustBePositive}
                ncz (1,1) {mustBeInteger, mustBePositive}
                ncoef (1,1) {mustBeInteger, mustBePositive}
                dt (1,1) {mustBeReal, mustBePositive}
                spht (1,1) {mustBeReal, mustBeNonnegative}
                con (1,1) {mustBeReal, mustBeNonnegative}
                heat_src (1,1) {mustBeReal}
                dens (:,:,:) {mustBeReal, mustBeFinite}
                t (:,:,:) {mustBeReal, mustBeFinite}
                uf (:,:,:) {mustBeReal, mustBeFinite}
                vf (:,:,:) {mustBeReal, mustBeFinite}
                wf (:,:,:) {mustBeReal, mustBeFinite}
                coef_tensor (:,:,:,:) {mustBeReal, mustBeFinite}
            end

            % 获取对流项离散格式
            scheme = case_obj.scheme;
            % 系数索引
            id_aP = case_obj.id_aP;
            id_aE = case_obj.id_aE;
            id_aW = case_obj.id_aW;
            id_aN = case_obj.id_aN;
            id_aS = case_obj.id_aS;
            if dim == 3
                id_aT = case_obj.id_aT;
                id_aB = case_obj.id_aB;
            end
            id_bsrc = case_obj.id_bsrc;
            idt = 1.0 / dt;

            % 获取网格几何信息
            dx = case_obj.dx;
            dy = case_obj.dy;
            dz = case_obj.dz;
            [area_x, area_y, area_z, vol] = solver.cal_area_vol(dx, dy, dz);

            idx = 1.0 / dx;
            idy = 1.0 / dy;
            idz = 1.0 / dz;

            % 预分配系数数组
            aE = zeros(ncx, ncy, ncz);
            aW = zeros(ncx, ncy, ncz);
            aN = zeros(ncx, ncy, ncz);
            aS = zeros(ncx, ncy, ncz);
            if dim == 3
                aT = zeros(ncx, ncy, ncz);
                aB = zeros(ncx, ncy, ncz);
            end

            % 源项
            sC = heat_src;
            sP = 0;

            % 计算密度插值
            % 东方向
            rho_e = zeros(ncx, ncy, ncz);
            rho_e(1:end-1,:,:) = 0.5 * (dens(1:end-1,:,:) + dens(2:end,:,:));
            rho_e(end,:,:) = dens(end,:,:);

            % 西方向
            rho_w = zeros(ncx, ncy, ncz);
            rho_w(2:end,:,:) = 0.5 * (dens(1:end-1,:,:) + dens(2:end,:,:));
            rho_w(1,:,:) = dens(1,:,:);

            % 北方向
            rho_n = zeros(ncx, ncy, ncz);
            rho_n(:,1:end-1,:) = 0.5 * (dens(:,1:end-1,:) + dens(:,2:end,:));
            rho_n(:,end,:) = dens(:,end,:);

            % 南方向
            rho_s = zeros(ncx, ncy, ncz);
            rho_s(:,2:end,:) = 0.5 * (dens(:,1:end-1,:) + dens(:,2:end,:));
            rho_s(:,1,:) = dens(:,1,:);

            % 3D情况
            if dim == 3
                % 顶面
                rho_t = zeros(ncx, ncy, ncz);
                rho_t(:,:,1:end-1) = 0.5 * (dens(:,:,1:end-1) + dens(:,:,2:end));
                rho_t(:,:,end) = dens(:,:,end);

                % 底面
                rho_b = zeros(ncx, ncy, ncz);
                rho_b(:,:,2:end) = 0.5 * (dens(:,:,1:end-1) + dens(:,:,2:end));
                rho_b(:,:,1) = dens(:,:,1);
            end

            % 计算系数
            mul = con; mur = con;

    % 东方向系数
    ul = uf(2:end,:,:); ur = uf(2:end,:,:);
    aE(:,:,:) = solver.a_nb(area_x, idx, ul, ur, mul, mur, rho_e(:,:,:), -1.0, scheme);
    
    % 西方向系数
    ul = uf(1:end-1,:,:); ur = uf(1:end-1,:,:);
    aW(:,:,:) = solver.a_nb(area_x, idx, ul, ur, mul, mur, rho_w(:,:,:), 1.0, scheme);
    
    % 北方向系数
    vl = vf(:,2:end,:); vr = vf(:,2:end,:);
    aN(:,:,:) = solver.a_nb(area_y, idy, vl, vr, mul, mur, rho_n(:,:,:), -1.0, scheme);
    
    % 南方向系数  
    vl = vf(:,1:end-1,:); vr = vf(:,1:end-1,:);
    aS(:,:,:) = solver.a_nb(area_y, idy, vl, vr, mul, mur, rho_s(:,:,:), 1.0, scheme);
    
    if dim == 3
        % 顶面系数
        wl = wf(:,:,2:end); wr = wf(:,:,2:end);
        aT(:,:,:) = solver.a_nb(area_z, idz, wl, wr, mul, mur, rho_t(:,:,:), -1.0, scheme);
        
        % 底面系数
        wl = wf(:,:,1:end-1); wr = wf(:,:,1:end-1);
        aB(:,:,:) = solver.a_nb(area_z, idz, wl, wr, mul, mur, rho_b(:,:,:), 1.0, scheme);
    end

            % 计算aP0和总aP
            aP0 = dens .* spht .* vol .* idt;
            aP = aE + aW + aN + aS;
            if dim == 3
                aP = aP + aT + aB;
            end
            aP = aP + aP0 - sP*vol;

            bsrc = sC .* vol + aP0 .* t;

            % 存储系数
            coef_tensor(:,:,:,id_aP) = aP;
            coef_tensor(:,:,:,id_aE) = -aE;
            coef_tensor(:,:,:,id_aW) = -aW;
            coef_tensor(:,:,:,id_aN) = -aN;
            coef_tensor(:,:,:,id_aS) = -aS;
            if dim == 3
                coef_tensor(:,:,:,id_aT) = -aT;
                coef_tensor(:,:,:,id_aB) = -aB;
            end
            coef_tensor(:,:,:,id_bsrc) = bsrc;
        end

        function coef_tensor = conduction_coef_bound(case_obj, boundary_obj, dim, ncx, ncy, ncz, ...
        ncoef, dt, con, dens, t, coef_tensor)
            % 处理边界条件
            %
            % 输入参数:
            %   case_obj - 案例对象
            %   boundary_obj - 流体边界对象
            %   dim - 维度
            %   ncx, ncy, ncz - 各方向单元数
            %   ncoef - 系数个数
            %   dt - 时间步长
            %   con - 传导系数
            %   dens - 密度场
            %   t - 温度场
            %   coef_tensor - 系数张量

            arguments
                case_obj (1,1) structured_mesh
                boundary_obj (1,1) boundary
                dim (1,1) {mustBeMember(dim, [2, 3])}
                ncx (1,1) {mustBeInteger, mustBePositive}
                ncy (1,1) {mustBeInteger, mustBePositive}
                ncz (1,1) {mustBeInteger, mustBePositive}
                ncoef (1,1) {mustBeInteger, mustBePositive}
                dt (1,1) {mustBeReal, mustBePositive}
                con (1,1) {mustBeReal, mustBeNonnegative}
                dens (:,:,:) {mustBeReal, mustBeFinite}
                t (:,:,:) {mustBeReal, mustBeFinite}
                coef_tensor (:,:,:,:) {mustBeReal, mustBeFinite}
            end

            % 系数索引
            id_aP = case_obj.id_aP;
            id_aE = case_obj.id_aE;
            id_aW = case_obj.id_aW;
            id_aN = case_obj.id_aN;
            id_aS = case_obj.id_aS;
            if dim == 3
                id_aT = case_obj.id_aT;
                id_aB = case_obj.id_aB;
            end
            id_bsrc = case_obj.id_bsrc;

            % 网格坐标
            x = case_obj.x;
            xc = case_obj.xc;
            y = case_obj.y;
            yc = case_obj.yc;
            if dim == 3
                z = case_obj.z;
                zc = case_obj.zc;
            end

            % 边界条件标识
            bcid = boundary_obj.bcid;

            % 面标识符
            fid_e = boundary_obj.fid_e;
            fid_w = boundary_obj.fid_w;
            fid_n = boundary_obj.fid_n;
            fid_s = boundary_obj.fid_s;
            if dim == 3
                fid_t = boundary_obj.fid_t;
                fid_b = boundary_obj.fid_b;
            end

            % 边界条件类型
            bcs = boundary_obj.bcs_fluid;
            bcs_temp = boundary_obj.bcs_temp;

            % bc_none   = boundary_obj.bc_none;
            bc_wall   = boundary_obj.bc_wall;
            bc_inlet  = boundary_obj.bc_inlet;
            bc_outlet = boundary_obj.bc_outlet;

            temp_bc_constant = boundary_obj.temp_bc_constant;
            temp_bc_heatflux = boundary_obj.temp_bc_heatflux;
            % % 对边界单元应用边界条件
            % 针对格式化网格的特点分别对单元的各个边界单独处理
            % 1. 西边界 (i=1)
            i = 1;
            for k = 1:ncz
                for j = 1:ncy
    
                    bcid_w = bcid(i,j,k,fid_w);
                    if bcid_w ~= 0
                        bc_w = bcs{bcid_w}.type;
                        temp_type_w = bcs_temp{bcid_w}.temp_type;
                        bc_temp_w = bcs_temp{bcid_w}.t;
                        bc_flux_w = bcs_temp{bcid_w}.heat_flux;
    
                        if bc_w == bc_wall
                            if temp_type_w == temp_bc_constant
                                coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aW);
                                coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aW) * bc_temp_w;
                                coef_tensor(i,j,k,id_aW) = 0.0;
                            elseif temp_type_w == temp_bc_heatflux
                                coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aW);
                                coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aW) * bc_flux_w * 2.0 * (x(i) - xc(i)) / con;
                                coef_tensor(i,j,k,id_aW) = 0.0;
                            end
                        end
                    end
                end
            end
    
            % 2. 东边界 (i=ncx)
            i = ncx;
            for k = 1:ncz
                for j = 1:ncy
    
                    bcid_e = bcid(i,j,k,fid_e);
                    if bcid_e ~= 0
                        bc_e = bcs{bcid_e}.type;
                        temp_type_e = bcs_temp{bcid_e}.temp_type;
                        bc_temp_e = bcs_temp{bcid_e}.t;
                        bc_flux_e = bcs_temp{bcid_e}.heat_flux;
    
                        if bc_e == bc_wall
                            if temp_type_e == temp_bc_constant
                                coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aE);
                                coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aE) * bc_temp_e;
                                coef_tensor(i,j,k,id_aE) = 0.0;
                            elseif temp_type_e == temp_bc_heatflux
                                coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aE);
                                coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aE) * bc_flux_e * 2.0 * (x(i+1) - xc(i)) / con;
                                coef_tensor(i,j,k,id_aE) = 0.0;
                            end
                        end
                    end
                end
            end
    
            % 3. 南边界 (j=1)
            j = 1;
            for k = 1:ncz
                for i = 1:ncx
    
                    bcid_s = bcid(i,j,k,fid_s);
                    if bcid_s ~= 0
                        bc_s = bcs{bcid_s}.type;
                        temp_type_s = bcs_temp{bcid_s}.temp_type;
                        bc_temp_s = bcs_temp{bcid_s}.t;
                        bc_flux_s = bcs_temp{bcid_s}.heat_flux;
    
                        if bc_s == bc_wall
                            if temp_type_s == temp_bc_constant
                                coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aS);
                                coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aS) * bc_temp_s;
                                coef_tensor(i,j,k,id_aS) = 0.0;
                            elseif temp_type_s == temp_bc_heatflux
                                coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aS);
                                coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aS) * bc_flux_s * 2.0 * (y(j) - yc(j)) / con;
                                coef_tensor(i,j,k,id_aS) = 0.0;
                            end
                        end
                    end
                end
            end
    
            % 4. 北边界 (j=ncy)
            j = ncy;
            for k = 1:ncz
                for i = 1:ncx
    
                    bcid_n = bcid(i,j,k,fid_n);
                    if bcid_n ~= 0
                        bc_n = bcs{bcid_n}.type;
                        temp_type_n = bcs_temp{bcid_n}.temp_type;
                        bc_temp_n = bcs_temp{bcid_n}.t;
                        bc_flux_n = bcs_temp{bcid_n}.heat_flux;
    
                        if bc_n == bc_wall
                            if temp_type_n == temp_bc_constant
                                coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aN);
                                coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aN) * bc_temp_n;
                                coef_tensor(i,j,k,id_aN) = 0.0;
                            elseif temp_type_n == temp_bc_heatflux
                                coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aN);
                                coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aN) * bc_flux_n * 2.0 * (y(j+1) - yc(j)) / con;
                                coef_tensor(i,j,k,id_aN) = 0.0;
                            end
                        end
                    end
                end
            end
    
            if dim == 3
            % 5. 底边界 (k=1) - 仅3D情况
                k = 1;
                for j = 1:ncy
                    for i = 1:ncx
    
                        bcid_b = bcid(i,j,k,fid_b);
                        if bcid_b ~= 0
                            bc_b = bcs{bcid_b}.type;
                            temp_type_b = bcs_temp{bcid_b}.temp_type;
                            bc_temp_b = bcs_temp{bcid_b}.t;
                            bc_flux_b = bcs_temp{bcid_b}.heat_flux;
    
                            if bc_b == bc_wall
                                if temp_type_b == temp_bc_constant
                                    coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aB);
                                    coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aB) * bc_temp_b;
                                    coef_tensor(i,j,k,id_aB) = 0.0;
                                elseif temp_type_b == temp_bc_heatflux
                                    coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aB);
                                    coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aB) * bc_flux_b * 2.0 * (z(k) - zc(k)) / con;
                                    coef_tensor(i,j,k,id_aB) = 0.0;
                                end
                            end
                        end
                    end
                end
    
            % 6. 顶边界 (k=ncz) - 仅3D情况
                k = ncz;
                for j = 1:ncy
                    for i = 1:ncx
    
                        bcid_t = bcid(i,j,k,fid_t);
                        if bcid_t ~= 0
                            bc_t = bcs{bcid_t}.type;
                            temp_type_t = bcs_temp{bcid_t}.temp_type;
                            bc_temp_t = bcs_temp{bcid_t}.t;
                            bc_flux_t = bcs_temp{bcid_t}.heat_flux;
    
                            if bc_t == bc_wall
                                if temp_type_t == temp_bc_constant
                                    coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aT);
                                    coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aT) * bc_temp_t;
                                    coef_tensor(i,j,k,id_aT) = 0.0;
                                elseif temp_type_t == temp_bc_heatflux
                                    coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aT);
                                    coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aT) * bc_flux_t * 2.0 * (z(k+1) - zc(k)) / con;
                                    coef_tensor(i,j,k,id_aT) = 0.0;
                                end
                            end
                        end
                    end
                end
            end
        end

        % function coef_tensor = conduction_coef_bound(case_obj, boundary_obj, dim, ncx, ncy, ncz, ...
        % ncoef, dt, con, dens, t, coef_tensor)
        %     % 处理边界条件
        %     %
        %     % 输入参数:
        %     %   case_obj - 案例对象
        %     %   boundary_obj - 流体边界对象
        %     %   dim - 维度
        %     %   ncx, ncy, ncz - 各方向单元数
        %     %   ncoef - 系数个数
        %     %   dt - 时间步长
        %     %   con - 传导系数
        %     %   dens - 密度场
        %     %   t - 温度场
        %     %   coef_tensor - 系数张量
        % 
        %     arguments
        %         case_obj (1,1) structured_mesh
        %         boundary_obj (1,1) boundary
        %         dim (1,1) {mustBeMember(dim, [2, 3])}
        %         ncx (1,1) {mustBeInteger, mustBePositive}
        %         ncy (1,1) {mustBeInteger, mustBePositive}
        %         ncz (1,1) {mustBeInteger, mustBePositive}
        %         ncoef (1,1) {mustBeInteger, mustBePositive}
        %         dt (1,1) {mustBeReal, mustBePositive}
        %         con (1,1) {mustBeReal, mustBeNonnegative}
        %         dens (:,:,:) {mustBeReal, mustBeFinite}
        %         t (:,:,:) {mustBeReal, mustBeFinite}
        %         coef_tensor (:,:,:,:) {mustBeReal, mustBeFinite}
        %     end
        % 
        %     % 系数索引
        %     id_aP = case_obj.id_aP;
        %     id_aE = case_obj.id_aE;
        %     id_aW = case_obj.id_aW;
        %     id_aN = case_obj.id_aN;
        %     id_aS = case_obj.id_aS;
        %     if dim == 3
        %         id_aT = case_obj.id_aT;
        %         id_aB = case_obj.id_aB;
        %     end
        %     id_bsrc = case_obj.id_bsrc;
        % 
        %     % 网格坐标
        %     x = case_obj.x;
        %     xc = case_obj.xc;
        %     y = case_obj.y;
        %     yc = case_obj.yc;
        %     if dim == 3
        %         z = case_obj.z;
        %         zc = case_obj.zc;
        %     end
        % 
        %     % 边界条件标识
        %     bcid = boundary_obj.bcid;
        % 
        %     % 面标识符
        %     fid_e = boundary_obj.fid_e;
        %     fid_w = boundary_obj.fid_w;
        %     fid_n = boundary_obj.fid_n;
        %     fid_s = boundary_obj.fid_s;
        %     if dim == 3
        %         fid_t = boundary_obj.fid_t;
        %         fid_b = boundary_obj.fid_b;
        %     end
        % 
        %     % 边界条件类型
        %     bcs = boundary_obj.bcs_fluid;
        %     bcs_temp = boundary_obj.bcs_temp;
        % 
        %     % bc_none   = boundary_obj.bc_none;
        %     bc_wall   = boundary_obj.bc_wall;
        %     bc_inlet  = boundary_obj.bc_inlet;
        %     bc_outlet = boundary_obj.bc_outlet;
        % 
        %     temp_bc_constant = boundary_obj.temp_bc_constant;
        %     temp_bc_heatflux = boundary_obj.temp_bc_heatflux;
        % 
        %     % 对每个单元应用边界条件
        %     for k = 1:ncz
        %         for j = 1:ncy
        %             for i = 1:ncx
        %                 % 获取边界标识 - 只在非零时才处理
        %                 bcid_e = bcid(i,j,k,fid_e);
        %                 bcid_w = bcid(i,j,k,fid_w);
        %                 bcid_n = bcid(i,j,k,fid_n);
        %                 bcid_s = bcid(i,j,k,fid_s);
        % 
        %                 % 东方向边界处理 - 只在有边界时才获取边界条件类型
        %                 if bcid_e ~= 0
        %                     bc_e = bcs{bcid_e}.type;
        %                     temp_type_e = bcs_temp{bcid_e}.temp_type;
        %                     bc_temp_e = bcs_temp{bcid_e}.t;
        %                     bc_flux_e = bcs_temp{bcid_e}.heat_flux;
        % 
        %                     if bc_e == bc_inlet
        %                         coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aE);
        %                         coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aE) * bc_temp_e;
        %                         coef_tensor(i,j,k,id_aE) = 0.0;
        %                     elseif bc_e == bc_outlet
        %                         coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aE);
        %                         coef_tensor(i,j,k,id_aE) = 0.0;
        %                     elseif bc_e == bc_wall
        %                         if temp_type_e == temp_bc_constant
        %                             coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aE);
        %                             coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aE) * bc_temp_e;
        %                             coef_tensor(i,j,k,id_aE) = 0.0;
        %                         elseif temp_type_e == temp_bc_heatflux
        %                             coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aE);
        %                             coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aE) * bc_flux_e * 2.0 * (x(i+1) - xc(i)) / con;
        %                             coef_tensor(i,j,k,id_aE) = 0.0;
        %                         end
        %                     end
        %                 end
        % 
        %                 % 西方向边界处理 - 只在有边界时才获取边界条件类型
        %                 if bcid_w ~= 0
        %                     bc_w = bcs{bcid_w}.type;
        %                     temp_type_w = bcs_temp{bcid_w}.temp_type;
        %                     bc_temp_w = bcs_temp{bcid_w}.t;
        %                     bc_flux_w = bcs_temp{bcid_w}.heat_flux;
        % 
        %                     if bc_w == bc_inlet
        %                         coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aW);
        %                         coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aW) * bc_temp_w;
        %                         coef_tensor(i,j,k,id_aW) = 0.0;
        %                     elseif bc_w == bc_outlet
        %                         coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aW);
        %                         coef_tensor(i,j,k,id_aW) = 0.0;
        %                     elseif bc_w == bc_wall
        %                         if temp_type_w == temp_bc_constant
        %                             coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aW);
        %                             coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aW) * bc_temp_w;
        %                             coef_tensor(i,j,k,id_aW) = 0.0;
        %                         elseif temp_type_w == temp_bc_heatflux
        %                             coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aW);
        %                             coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aW) * bc_flux_w * 2.0 * (x(i) - xc(i)) / con;
        %                             coef_tensor(i,j,k,id_aW) = 0.0;
        %                         end
        %                     end
        %                 end
        % 
        %                 % 北方向边界处理 - 只在有边界时才获取边界条件类型
        %                 if bcid_n ~= 0
        %                     bc_n = bcs{bcid_n}.type;
        %                     temp_type_n = bcs_temp{bcid_n}.temp_type;
        %                     bc_temp_n = bcs_temp{bcid_n}.t;
        %                     bc_flux_n = bcs_temp{bcid_n}.heat_flux;
        % 
        %                     if bc_n == bc_inlet
        %                         coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aN);
        %                         coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aN) * bc_temp_n;
        %                         coef_tensor(i,j,k,id_aN) = 0.0;
        %                     elseif bc_n == bc_outlet
        %                         coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aN);
        %                         coef_tensor(i,j,k,id_aN) = 0.0;
        %                     elseif bc_n == bc_wall
        %                         if temp_type_n == temp_bc_constant
        %                             coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aN);
        %                             coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aN) * bc_temp_n;
        %                             coef_tensor(i,j,k,id_aN) = 0.0;
        %                         elseif temp_type_n == temp_bc_heatflux
        %                             coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aN);
        %                             coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aN) * bc_flux_n * 2.0 * (y(j+1) - yc(j)) / con;
        %                             coef_tensor(i,j,k,id_aN) = 0.0;
        %                         end
        %                     end
        %                 end
        % 
        %                 % 南方向边界处理 - 只在有边界时才获取边界条件类型
        %                 if bcid_s ~= 0
        %                     bc_s = bcs{bcid_s}.type;
        %                     temp_type_s = bcs_temp{bcid_s}.temp_type;
        %                     bc_temp_s = bcs_temp{bcid_s}.t;
        %                     bc_flux_s = bcs_temp{bcid_s}.heat_flux;
        % 
        %                     if bc_s == bc_inlet
        %                         coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aS);
        %                         coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aS) * bc_temp_s;
        %                         coef_tensor(i,j,k,id_aS) = 0.0;
        %                     elseif bc_s == bc_outlet
        %                         coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aS);
        %                         coef_tensor(i,j,k,id_aS) = 0.0;
        %                     elseif bc_s == bc_wall
        %                         if temp_type_s == temp_bc_constant
        %                             coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) - coef_tensor(i,j,k,id_aS);
        %                             coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) - 2.0 * coef_tensor(i,j,k,id_aS) * bc_temp_s;
        %                             coef_tensor(i,j,k,id_aS) = 0.0;
        %                         elseif temp_type_s == temp_bc_heatflux
        %                             coef_tensor(i,j,k,id_aP) = coef_tensor(i,j,k,id_aP) + coef_tensor(i,j,k,id_aS);
        %                             coef_tensor(i,j,k,id_bsrc) = coef_tensor(i,j,k,id_bsrc) + coef_tensor(i,j,k,id_aS) * bc_flux_s * 2.0 * (y(j) - yc(j)) / con;
        %                             coef_tensor(i,j,k,id_aS) = 0.0;
        %                         end
        %                     end
        %                 end
        %             end
        %         end
        %     end

        


    end
end

