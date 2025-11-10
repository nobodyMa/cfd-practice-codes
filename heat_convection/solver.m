classdef solver
    methods (Static)
        function [area_x, area_y, area_z, vol] = cal_area_vol(dx, dy, dz)
            % 计算面积和体积
            area_x = dy * dz;
            area_y = dx * dz;
            area_z = dx * dy;
            vol = dx * area_x;
        end

        function a = a_nb(area, idx, ul, ur, gl, gr, rho, sign_f, scheme)
            % 计算邻居系数
            arguments
                area
                idx
                ul
                ur
                gl
                gr
                rho
                sign_f
                scheme {mustBeMember(scheme, {'upwind', 'centraldifference'})} = 'upwind'
            end
            f = 0.5 .* rho .* (ul + ur);
            d = 2.0 .* gl .* gr ./ (gl + gr + 1.e-12) .* idx;
            switch scheme
                case 'upwind'
                    % 迎风格式
                    a = area .* (d + max(0, sign_f .* f));
                case 'centraldifference'
                    % 中心差分格式 - 对流项使用中心差分
                    a = area .* (d .* (1.0 - 0.5 .* abs(f/d)) + max(0, sign_f*f));
                otherwise
                    error('不支持的格式: %s', scheme);
            end
        end

        function [t, rel_norm] = scalar_pj(case_obj, post, dim, it_nl, niter, relax, ncx, ncy, ncz, ncoef, ct, t, initzero, res)
            % 雅可比迭代法 - 向量化运算
            % case_obj - 案例对象
            % post - 后处理对象
            % dim - 维度
            % it_nl - 非线性迭代次数
            % niter - 线性求解器迭代次数
            % relax - 松弛因子
            % ncx, ncy, ncz - 各方向单元数
            % ncoef - 系数个数
            % ct - 系数张量
            % t - 温度场
            % initzero - 是否初始化为零
            % res - 收敛容差

            % 预计算系数矩阵
            aP = ct(:,:,:,case_obj.id_aP);
            aE = ct(:,:,:,case_obj.id_aE);
            aW = ct(:,:,:,case_obj.id_aW);
            aN = ct(:,:,:,case_obj.id_aN);
            aS = ct(:,:,:,case_obj.id_aS);
            bsrc = ct(:,:,:,case_obj.id_bsrc);
            
            if dim == 3
                aT = ct(:,:,:,case_obj.id_aT);
                aB = ct(:,:,:,case_obj.id_aB);
            end
            
            max_norm = -1e20;
            
            for it = 1:niter
                t0 = t;
                
                % 使用矩阵运算计算邻居贡献
                t_temp = bsrc;
                
                % 东方向邻居贡献
                t_temp(1:end-1,:,:) = t_temp(1:end-1,:,:) - aE(1:end-1,:,:) .* t0(2:end,:,:);
                % 西方向邻居贡献
                t_temp(2:end,:,:) = t_temp(2:end,:,:) - aW(2:end,:,:) .* t0(1:end-1,:,:);
                % 北方向邻居贡献
                t_temp(:,1:end-1,:) = t_temp(:,1:end-1,:) - aN(:,1:end-1,:) .* t0(:,2:end,:);
                % 南方向邻居贡献
                t_temp(:,2:end,:) = t_temp(:,2:end,:) - aS(:,2:end,:) .* t0(:,1:end-1,:);
                
                if dim == 3
                    % 顶方向邻居贡献
                    t_temp(:,:,1:end-1) = t_temp(:,:,1:end-1) - aT(:,:,1:end-1) .* t0(:,:,2:end);
                    % 底方向邻居贡献
                    t_temp(:,:,2:end) = t_temp(:,:,2:end) - aB(:,:,2:end) .* t0(:,:,1:end-1);
                end
                
                % 计算新温度
                t_new = t_temp ./ aP;
                
                % 应用松弛因子
                du = relax * (t_new - t0);
                t = t0 + du;
                
                % 更新最大范数
                if it == 1
                    max_norm = max(max_norm, norm(du(:)) / sqrt(numel(t)));
                end
            end
            
            % 计算范数
            norm2 = norm(du(:)) / sqrt(numel(t));
            rel_norm = norm2 / max_norm;
            case_obj.total_linsol_iters = case_obj.total_linsol_iters + it;
            
            % 输出收敛信息
            if it_nl == 1 || it_nl == case_obj.nsteps || mod(it_nl, post.res_freq) == 0
                fprintf('it_nl, it, tot_it, norm2, max, rel_norm, : %d, %d, %d, %.6e, %.6e, %.6e, %.6e\n', ...
                    it_nl, it, case_obj.total_linsol_iters, norm2, max_norm, rel_norm);
                
                % 写入线性求解器结果文件
                linsol_fid = fopen(post.linsol_fname, 'a');
                fprintf(linsol_fid, '%d %d %d %.6e %.6e %.6e\n', ...
                    it_nl, it, case_obj.total_linsol_iters, norm2, max_norm, rel_norm);
                fclose(linsol_fid);
            end
        end

    end
end