classdef post_processing
    properties
        % 输出频率
        res_freq = 1000
        out_freq = 2000
                
        % 输出文件名
        linsol_fname    = "lin.res"
        nonlinsol_fname = "nonlin.res"
        vtk_fname_temp  = "post_temp.vtk"

        % 输出文件句柄
        linsol_fid    = 101
        nonlinsol_fid = 102
        vtk_fid       = 103
        
    end
    
    methods
        function obj = post_processing()
            % 构造函数
        end
        
        function obj = Set_res_out_freq(obj, res_freq, out_freq)
            % 设置输出频率
            obj.res_freq = res_freq;
            obj.out_freq = out_freq;
        end
        
        function plot_temperature_2D(obj, mesh)
            % 绘制2D温度场
            
            if mesh.dim ~= 2
                error('此方法仅适用于2D网格');
            end
            
            figure;
            [X, Y] = meshgrid(mesh.xc, mesh.yc);
            
            % 转置温度矩阵以正确显示
            T_plot = mesh.t';
            
            contourf(X, Y, T_plot, 20, 'LineStyle', 'none');
            colorbar;
            xlabel('X');
            ylabel('Y');
            title('2D温度场分布');
            axis equal;
        end
        
        function plot_temperature_line(obj, mesh)
            % 绘制沿y线的温度分布
            
            % 读取温度数据文件
            data = load('temp_yline.dat');
            y = data(:, 1);
            t1 = data(:, 2);  % x = 0.0833 处的温度
            t2 = data(:, 3);  % x = 0.5 处的温度
            
            figure;
            plot(y, t1, 'r-', 'LineWidth', 2, 'DisplayName', 'x = 0.0833');
            hold on;
            plot(y, t2, 'b-', 'LineWidth', 2, 'DisplayName', 'x = 0.5');
            xlabel('Y坐标');
            ylabel('温度');
            title('沿Y线的温度分布');
            legend('Location', 'best');
            grid on;
        end
        
        function output_temperature_profile(obj, mesh, position_type, position_value, output_filename)
            % 输出任意位置的线或面的温度分布
            % position_type: 'x', 'y', 或 'z'
            % position_value: 位置坐标值
            % output_filename: 输出文件名
            
            arguments
                obj % 留给未来使用索引传递选择单元
                mesh (1,1) structured_mesh
                position_type {mustBeMember(position_type, {'x', 'y', 'z'})}
                position_value (1,1) {mustBeReal}
                output_filename {mustBeText} = 'temp_profile.dat'
            end
            
            % 获取网格坐标
            xc = mesh.xc;
            yc = mesh.yc;
            zc = mesh.zc;
            t = mesh.t;
            
            % 根据维度处理
            if mesh.dim == 2
                % 2D情况 - 输出线数据
                switch position_type
                    case 'x'
                        % 在指定x位置输出沿y方向的温度分布
                        [~, idx] = min(abs(xc - position_value));
                        y_coords = yc;
                        temp_values = squeeze(t(idx, :, 1));
                        
                        % 写入文件
                        fid = fopen(output_filename, 'w');
                        fprintf(fid, '# Y坐标, 温度 (x = %.6f)\n', xc(idx));
                        for j = 1:length(y_coords)
                            fprintf(fid, '%.6f %.6f\n', y_coords(j), temp_values(j));
                        end
                        fclose(fid);
                        
                        % 绘制图形
                        figure;
                        plot(y_coords, temp_values, 'b-', 'LineWidth', 2);
                        xlabel('Y坐标');
                        ylabel('温度');
                        title(sprintf('沿Y方向的温度分布 (x = %.4f)', xc(idx)));
                        grid on;
                        
                    case 'y'
                        % 在指定y位置输出沿x方向的温度分布
                        [~, idx] = min(abs(yc - position_value));
                        x_coords = xc;
                        temp_values = squeeze(t(:, idx, 1));
                        
                        % 写入文件
                        fid = fopen(output_filename, 'w');
                        fprintf(fid, '# X坐标, 温度 (y = %.6f)\n', yc(idx));
                        for i = 1:length(x_coords)
                            fprintf(fid, '%.6f %.6f\n', x_coords(i), temp_values(i));
                        end
                        fclose(fid);
                        
                        % 绘制图形
                        figure;
                        plot(x_coords, temp_values, 'r-', 'LineWidth', 2);
                        xlabel('X坐标');
                        ylabel('温度');
                        title(sprintf('沿X方向的温度分布 (y = %.4f)', yc(idx)));
                        grid on;
                        
                    otherwise
                        error('2D情况下不支持z方向');
                end
                
            else
                % 3D情况 - 输出面数据
                switch position_type
                    case 'x'
                        % 在指定x位置输出yz平面的温度分布
                        [~, idx] = min(abs(xc - position_value));
                        [Y, Z] = meshgrid(yc, zc);
                        temp_slice = squeeze(t(idx, :, :));
                        
                        % 写入文件
                        fid = fopen(output_filename, 'w');
                        fprintf(fid, '# Y坐标, Z坐标, 温度 (x = %.6f)\n', xc(idx));
                        for k = 1:length(zc)
                            for j = 1:length(yc)
                                fprintf(fid, '%.6f %.6f %.6f\n', yc(j), zc(k), temp_slice(j, k));
                            end
                            fprintf(fid, '\n');  % 空行分隔不同z值的数据
                        end
                        fclose(fid);
                        
                        % 绘制图形
                        figure;
                        contourf(Y, Z, temp_slice', 20, 'LineStyle', 'none');
                        colorbar;
                        xlabel('Y坐标');
                        ylabel('Z坐标');
                        title(sprintf('YZ平面温度分布 (x = %.4f)', xc(idx)));
                        axis equal;
                        
                    case 'y'
                        % 在指定y位置输出xz平面的温度分布
                        [~, idx] = min(abs(yc - position_value));
                        [X, Z] = meshgrid(xc, zc);
                        temp_slice = squeeze(t(:, idx, :));
                        
                        % 写入文件
                        fid = fopen(output_filename, 'w');
                        fprintf(fid, '# X坐标, Z坐标, 温度 (y = %.6f)\n', yc(idx));
                        for k = 1:length(zc)
                            for i = 1:length(xc)
                                fprintf(fid, '%.6f %.6f %.6f\n', xc(i), zc(k), temp_slice(i, k));
                            end
                            fprintf(fid, '\n');  % 空行分隔不同z值的数据
                        end
                        fclose(fid);
                        
                        % 绘制图形
                        figure;
                        contourf(X, Z, temp_slice', 20, 'LineStyle', 'none');
                        colorbar;
                        xlabel('X坐标');
                        ylabel('Z坐标');
                        title(sprintf('XZ平面温度分布 (y = %.4f)', yc(idx)));
                        axis equal;
                        
                    case 'z'
                        % 在指定z位置输出xy平面的温度分布
                        [~, idx] = min(abs(zc - position_value));
                        [X, Y] = meshgrid(xc, yc);
                        temp_slice = squeeze(t(:, :, idx));
                        
                        % 写入文件
                        fid = fopen(output_filename, 'w');
                        fprintf(fid, '# X坐标, Y坐标, 温度 (z = %.6f)\n', zc(idx));
                        for j = 1:length(yc)
                            for i = 1:length(xc)
                                fprintf(fid, '%.6f %.6f %.6f\n', xc(i), yc(j), temp_slice(i, j));
                            end
                            fprintf(fid, '\n');  % 空行分隔不同y值的数据
                        end
                        fclose(fid);
                        
                        % 绘制图形
                        figure;
                        contourf(X, Y, temp_slice', 20, 'LineStyle', 'none');
                        colorbar;
                        xlabel('X坐标');
                        ylabel('Y坐标');
                        title(sprintf('XY平面温度分布 (z = %.4f)', zc(idx)));
                        axis equal;
                end
            end
            fprintf('温度分布已输出到文件: %s\n', output_filename);
        end

        function plot_temperature_3D_dense_slices(obj, mesh, num_slices)
            % 使用密集切片体绘制显示3D温度场
            % num_slices: 每个方向的切片数量
            
            if mesh.dim ~= 3
                error('此方法仅适用于3D网格');
            end
            
            if nargin < 3
                num_slices = ceil(mesh.ncx / 4) + 1;  % 默认每个方向20个切片
            end
            
            % 获取网格数据
            xc = mesh.xc;
            yc = mesh.yc;
            zc = mesh.zc;
            t = mesh.t;
            
            % 创建3D网格
            [X, Y, Z] = meshgrid(xc, yc, zc);
            T_3D = permute(t, [2, 1, 3]);
            
            figure;
            hold on;
            % 生成密集切片位置
            x_slices = linspace(min(xc), max(xc), num_slices);
            y_slices = linspace(min(yc), max(yc), num_slices);
            z_slices = linspace(min(zc), max(zc), num_slices);
            
            % 创建密集切片
            hx = slice(X, Y, Z, T_3D, x_slices, [], []);
            set(hx, 'EdgeColor', 'none', 'FaceAlpha', 0.7);
            colorbar;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            axis equal; view(3); grid on;

            hy = slice(X, Y, Z, T_3D, [], y_slices, []);
            set(hy, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

            hz = slice(X, Y, Z, T_3D, [], [], z_slices);
            set(hz, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

            title(sprintf('每个方向密集切片 (%d个)', num_slices));
            hold off;
            fprintf('密集切片绘制完成 - 每个方向%d个切片\n', num_slices);
        end

    end
end