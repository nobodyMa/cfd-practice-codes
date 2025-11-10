classdef boundary
    properties
        mesh
        % 面标识符
        fid_e = 2      % 东面
        fid_w = 1      % 西面  
        fid_n = 4      % 北面
        fid_s = 3      % 南面
        fid_b = 5      % 底面
        fid_t = 6      % 顶面
        
        % 流体边界条件类型
        bc_none   = 0  % 无边界
        bc_wall   = 1  % 壁面
        bc_inlet  = 2  % 进口
        bc_outlet = 3  % 出口

        % 温度边界条件类型
        temp_bc_none = 0  % 无边界
        temp_bc_constant = 1  % 迪利克雷边界（定温）
        temp_bc_heatflux = 2  % 纽曼边界（热流）
        
        % 边界标识符
        % bcid_none = 0  % 非边界
        bcid_xmin = 1  % x最小边界
        bcid_xmax = 2  % x最大边界
        bcid_ymin = 3  % y最小边界
        bcid_ymax = 4  % y最大边界
        bcid_zmin = 5  % z最小边界
        bcid_zmax = 6  % z最大边界
        num_bcs   = 0  % 边界类型总数
        
        % 边界条件数组
        bcs_fluid       % 流体边界条件
        bcs_temp  % 温度边界条件
        
        % 温度边界条件值
        temp_east   = 0.0
        temp_west   = 0.0  
        temp_north  = 0.0
        temp_south  = 0.0
        temp_top    = 0.0
        temp_bottom = 0.0
               
        % 边界标识数组
        bcid      % 单元面边界标识 [ncx, ncy, ncz, 6]
    end
    
    methods
        function obj = boundary(mesh)
            % 构造函数 - 初始化边界条件数组
            obj.num_bcs = 2 * mesh.dim;
            obj.bcs_fluid = cell(obj.num_bcs, 1);
            obj.bcs_temp = cell(obj.num_bcs, 1);
            
            for i = 1:obj.num_bcs
                obj.bcs_fluid{i} = struct('type', obj.bc_none);
                obj.bcs_temp{i} = struct('temp_type', obj.temp_bc_none, ...
                                        't', 0.0, ...
                                        'heat_flux', 0.0);
            end
            
        end
        
        function obj = create_faces_of_cells(obj, mesh)
            % 创建单元面边界标识
            arguments
                obj
                mesh (1,1) structured_mesh
            end
            
            % 从网格对象获取参数
            dim  = mesh.dim;
            ncx  = mesh.ncx;
            ncy  = mesh.ncy;
            ncz  = mesh.ncz;
            xmin = mesh.xmin;
            xmax = mesh.xmax;
            ymin = mesh.ymin;
            ymax = mesh.ymax;
            zmin = mesh.zmin;
            zmax = mesh.zmax;
            x    = mesh.x;
            y    = mesh.y;
            z    = mesh.z;
            
            % 初始化边界标识数组
            obj.bcid = zeros(ncx, ncy, ncz, obj.num_bcs, 'int32');
            
            % 容差
            eps = 1e-12;
            
            % 西边界 (xmin)
            west_mask = abs(x(1:ncx) - xmin) < eps;
            obj.bcid(west_mask, :, :, obj.fid_w) = obj.bcid_xmin;
            
            % 东边界 (xmax)
            east_mask = abs(x(2:ncx+1) - xmax) < eps;
            obj.bcid(east_mask, :, :, obj.fid_e) = obj.bcid_xmax;

            % 南边界 (ymin)
            south_mask = abs(y(1:ncy) - ymin) < eps;
            obj.bcid(:, south_mask, :, obj.fid_s) = obj.bcid_ymin;
            
            % 北边界 (ymax)
            north_mask = abs(y(2:ncy+1) - ymax) < eps;
            obj.bcid(:, north_mask, :, obj.fid_n) = obj.bcid_ymax;
            
            % 3D情况
            if dim == 3
                % 底边界 (zmin)
                bottom_mask = abs(z(1:ncz) - zmin) < eps;
                obj.bcid(:, :, bottom_mask, obj.fid_b) = obj.bcid_zmin;
                
                % 顶边界 (zmax)
                top_mask = abs(z(2:ncz+1) - zmax) < eps;
                obj.bcid(:, :, top_mask, obj.fid_t) = obj.bcid_zmax;
            end
        end
        
        function obj = create_boundary_fluid(obj, dim, bc_xmin, bc_xmax, bc_ymin, bc_ymax, bc_zmin, bc_zmax)
            % 创建边界物理条件数据
            % dim - 维度
            % bc_xmin, bc_xmax, bc_ymin, bc_ymax, bc_zmin, bc_zmax - 各边界条件类型
            
            arguments
                obj
                dim (1,1) {mustBeMember(dim, [2, 3])}
                bc_xmin {mustBeText} = 'wall'
                bc_xmax {mustBeText} = 'wall'
                bc_ymin {mustBeText} = 'wall'
                bc_ymax {mustBeText} = 'wall'
                bc_zmin {mustBeText} = 'wall'
                bc_zmax {mustBeText} = 'wall'
            end
                       
            % 设置各边界条件
            boundaries = {
                obj.bcid_xmin, bc_xmin;
                obj.bcid_xmax, bc_xmax;
                obj.bcid_ymin, bc_ymin;
                obj.bcid_ymax, bc_ymax
            };
            
            if dim == 3
                boundaries = [boundaries; {
                    obj.bcid_zmin, bc_zmin;
                    obj.bcid_zmax, bc_zmax
                }];
            end
            
            for i = 1:size(boundaries, 1)
                id = boundaries{i,1};
                bc_type = boundaries{i,2};
                
                switch bc_type
                    case 'wall'
                        obj.bcs_fluid{id}.type = obj.bc_wall;
                        
                    otherwise
                        error('不支持的边界条件类型: %s', bc_type);
                end
            end
        end
        
        function obj = create_boundary_temp(obj, dim, varargin)
            % 创建温度边界条件数据
            % 参数格式: 边界，{温度类型, 温度值/热流值}
            
            p = inputParser;
            addRequired(p, 'obj');
            addRequired(p, 'dim', @(x) ismember(x, [2,3]));
            addParameter(p, 'xmin', {'constant', 0.0});
            addParameter(p, 'xmax', {'constant', 0.0});
            addParameter(p, 'ymin', {'constant', 0.0});
            addParameter(p, 'ymax', {'constant', 0.0});
            addParameter(p, 'zmin', {'constant', 0.0});
            addParameter(p, 'zmax', {'constant', 0.0});
            parse(p, obj, dim, varargin{:});
                        
            % 设置各边界温度条件
            boundaries = {
                obj.bcid_xmin, p.Results.xmin;
                obj.bcid_xmax, p.Results.xmax;
                obj.bcid_ymin, p.Results.ymin;
                obj.bcid_ymax, p.Results.ymax
            };
            
            if dim == 3
                boundaries = [boundaries; {
                    obj.bcid_zmin, p.Results.zmin;
                    obj.bcid_zmax, p.Results.zmax
                }];
            end
            
            for i = 1:size(boundaries, 1)
                id = boundaries{i,1};
                bc_data = boundaries{i,2};

                switch bc_data{1}
                    case 'constant'
                        obj.bcs_temp{id}.temp_type = obj.temp_bc_constant;
                        obj.bcs_temp{id}.t = bc_data{2};
                        obj.bcs_temp{id}.heat_flux = 0.0;
                    case 'heatflux'
                        obj.bcs_temp{id}.temp_type = obj.temp_bc_heatflux;
                        obj.bcs_temp{id}.heat_flux = bc_data{2};
                        obj.bcs_temp{id}.t = 0.0;
                    otherwise
                        error('不支持的温度边界条件类型: %s', temp_type_str);
                end
            end

            % display_info(obj);
        end
        
        function display_info(obj)
            % 显示边界条件信息
            fprintf('=== 边界条件信息 ===\n');
            fprintf('边界标识符: \n');
            fprintf('  xmin: %d, xmax: %d\n', obj.bcid_xmin, obj.bcid_xmax);
            fprintf('  ymin: %d, ymax: %d\n', obj.bcid_ymin, obj.bcid_ymax);
            fprintf('  zmin: %d, zmax: %d\n', obj.bcid_zmin, obj.bcid_zmax);
            
            boundary_names = {'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'};
            fluid_type_names = {'无限制', '壁面', '进口', '出口'};
            temp_type_names = {'无限制', '定温', '热流'};
            
            fprintf('\n流体边界条件:\n');
            for i = 1:obj.num_bcs
                fprintf('  %s: %s\n', boundary_names{i}, fluid_type_names{obj.bcs_fluid{i}.type + 1});
            end
            
            fprintf('\n温度边界条件:\n');
            for i = 1:obj.num_bcs
                fprintf('  %s: %s', boundary_names{i}, temp_type_names{obj.bcs_temp{i}.temp_type + 1});
                if obj.bcs_temp{i}.temp_type == obj.temp_bc_constant
                    fprintf(' (T=%.3f)', obj.bcs_temp{i}.t);
                elseif obj.bcs_temp{i}.temp_type == obj.temp_bc_heatflux
                    fprintf(' (q=%.3f)', obj.bcs_temp{i}.heat_flux);
                end
                fprintf('\n');
            end
        end
    end
end