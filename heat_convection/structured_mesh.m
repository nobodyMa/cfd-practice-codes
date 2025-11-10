classdef structured_mesh < handle
    properties
        dim     % 空间维度 (2 或 3)
        ncx     % x方向单元数
        ncy     % y方向单元数  
        ncz     % z方向单元数
        nx      % x方向节点数
        ny      % y方向节点数
        nz      % z方向节点数
        
        % 域坐标范围
        xmin = 0.0
        xmax = 0.0
        ymin = 0.0
        ymax = 0.0
        zmin = 0.0
        zmax = 0.0
        
        % 网格坐标
        x = []
        y = []
        z = []
        
        % 单元中心坐标
        xc = []
        yc = []
        zc = []
        
        % 网格间距
        dx = 0.0
        dy = 0.0
        dz = 0.0
    
        % 场数据
        t = []      % 温度场
        t0 = []   % 旧温度场
        
        % 速度场
        uf = []     % x方向速度场
        vf = []     % y方向速度场  
        wf = []     % z方向速度场

        % 系数存储
        ncoef    % 系数个数
        id_aP    % aP系数索引
        id_aE    % aE系数索引  
        id_aW    % aW系数索引
        id_aN    % aN系数索引
        id_aS    % aS系数索引
        id_aT    % aT系数索引
        id_aB    % aB系数索引
        id_bsrc  % 源项系数索引，在matlab中id_aP=id_bsrc，但在其他语言中可能不相等
        
        % 温度方程系数
        coef     % 系数数组 [ncx, ncy, ncz, ncoef]
        
        % 仿真控制参数
        ieqn = 0                    % 方程类型
        eqn_conduction = 0          % 只求解温度场
        eqn_flow = 1                % 只求解流场
        eqn_conduction_flow = 2     % 同时求解流场和温度场
        nsteps = 1                  % 非线性迭代次数
        dt = 1.0                    % 时间步长
        stop_sim = false            % 仿真停止标志
        
        % 对流项离散格式
        scheme = 'upwind'

        % 求解器参数
        niter_t = 10                % 温度线性求解器迭代次数
        relax_t = 0.75              % 松弛因子
        res_t = 1.e-2               % 线性求解器收敛误差
        total_linsol_iters = 0      % 线性迭代总次数
        temp_tol = 1.e-6            % 非线性迭代收敛误差
        
        % 残差范数
        norm2_curr = 0.0             % 当前迭代残差2范数
        norm2_max = -1e20            % 最大残差2范数
        norm2_max_temp = -1e20       % 温度相关最大残差2范数
        norm2_temp = 0.0             % 温度残差2范数
    end


    methods        
        function obj = create_mesh(obj, dim, ncx, ncy, ncz, xmin, xmax, ymin, ymax, zmin, zmax)
            % 创建网格            
            arguments
                obj
                dim (1,1) {mustBeMember(dim, [2, 3])}
                ncx (1,1) {mustBePositive, mustBeInteger}
                ncy (1,1) {mustBePositive, mustBeInteger}
                ncz (1,1) {mustBePositive, mustBeInteger}
                xmin (1,1) {mustBeReal}
                xmax (1,1) {mustBeReal}
                ymin (1,1) {mustBeReal}
                ymax (1,1) {mustBeReal}
                % zmin和zmax有默认值
                zmin (1,1) {mustBeReal} = 0
                zmax (1,1) {mustBeReal} = 1
            end
            
            % 设置维度
            obj.dim = dim;
            
            % 设置单元数
            obj.ncx = ncx;
            obj.ncy = ncy;
            obj.ncz = ncz;
            
            if obj.dim == 2
                obj.ncz = 1;
                fprintf('调试：2D网格检测成功，ncz 已设置为 1\n');
            end
            
            % 计算节点数
            obj.nx = obj.ncx + 1;
            obj.ny = obj.ncy + 1;
            obj.nz = obj.ncz + 1;
            
            % 设置域坐标范围
            obj.xmin = xmin;
            obj.xmax = xmax;
            obj.ymin = ymin;
            obj.ymax = ymax;
            obj.zmin = zmin;
            obj.zmax = zmax;
            
            % 生成网格坐标
            obj = obj.generate_mesh_coordinates();
            % 显示网格信息
            % display_info(obj)
        end
        
        function obj = generate_mesh_coordinates(obj)
            % 网格坐标生成
            obj.dx = (obj.xmax - obj.xmin) / (obj.nx - 1);
            obj.x = obj.xmin:obj.dx:obj.xmax;
            
            obj.dy = (obj.ymax - obj.ymin) / (obj.ny - 1);
            obj.y = obj.ymin:obj.dy:obj.ymax;
            
            obj.dz = (obj.zmax - obj.zmin) / (obj.nz - 1);
            obj.z = obj.zmin:obj.dz:obj.zmax;
                    
            % 单元中心坐标生成
            obj.xc = zeros(1, obj.ncx);
            for i = 1:obj.ncx
                obj.xc(i) = (obj.x(i) + obj.x(i+1)) / 2.0;
            end
            
            obj.yc = zeros(1, obj.ncy);
            for i = 1:obj.ncy
                obj.yc(i) = (obj.y(i) + obj.y(i+1)) / 2.0;
            end
            
            obj.zc = zeros(1, obj.ncz);
            for i = 1:obj.ncz
                obj.zc(i) = (obj.z(i) + obj.z(i+1)) / 2.0;
            end
        end
                
        function display_info(obj)
            % 显示网格信息
            fprintf('=== 结构化网格信息 ===\n');
            fprintf('维度: %dD\n', obj.dim);
            fprintf('单元数: %d (x) × %d (y) × %d (z)\n', obj.ncx, obj.ncy, obj.ncz);
            fprintf('节点数: %d (x) × %d (y) × %d (z)\n', obj.nx, obj.ny, obj.nz);
            fprintf('网格间距: dx=%.6f, dy=%.6f, dz=%.6f\n', obj.dx, obj.dy, obj.dz);
            fprintf('计算域范围: x=[%.3f, %.3f], y=[%.3f, %.3f]', ...
                obj.xmin, obj.xmax, obj.ymin, obj.ymax);
            if obj.dim == 3
                fprintf(', z=[%.3f, %.3f]', obj.zmin, obj.zmax);
            end
            fprintf('\n');
        end

        function obj = create_field_data(obj)
            % 创建场数据（包含温度场和速度场）
            
            % 温度场
            obj.t = zeros(obj.ncx, obj.ncy, obj.ncz);
            obj.t0 = zeros(obj.ncx, obj.ncy, obj.ncz);
            
            % 网格面上的速度场
            obj.uf = zeros(obj.nx, obj.ncy, obj.ncz);   % x方向速度场
            obj.vf = zeros(obj.ncx, obj.ny, obj.ncz);   % y方向速度场  
            obj.wf = zeros(obj.ncx, obj.ncy, obj.nz);   % z方向速度场
        end
        
        function obj = set_initial_t(obj, T)
            % 设置初始温度场
            obj.t = T * ones(obj.ncx, obj.ncy, obj.ncz);
        end

        function obj = set_face_velocity_field(obj, u, v, w)
            % 直接设置面速度场
            % u: [nx, ncy, ncz] - x方向面速度
            % v: [ncx, ny, ncz] - y方向面速度
            % w: [ncx, ncy, nz] - z方向面速度
                       
            % 检查尺寸
            if ~isequal(size(u), [obj.nx, obj.ncy, obj.ncz])
                error('x方向速度场u的尺寸应为 [%d, %d, %d]', obj.nx, obj.ncy, obj.ncz);
            end
            if ~isequal(size(v), [obj.ncx, obj.ny, obj.ncz])
                error('y方向速度场v的尺寸应为 [%d, %d, %d]', obj.ncx, obj.ny, obj.ncz);
            end
            if obj.dim == 3 && ~isequal(size(w), [obj.ncx, obj.ncy, obj.nz])
                error('z方向速度场w的尺寸应为 [%d, %d, %d]', obj.ncx, obj.ncy, obj.nz);
            end
            
            obj.uf = u;
            obj.vf = v;
            obj.wf = w;
            
            fprintf('面速度场设置完成\n');
        end
        function obj = set_uniform_face_velocity(obj, u_val, v_val, w_val)
            % 设置均匀面速度场
            arguments
                obj
                u_val (1,1) {mustBeReal} = 0
                v_val (1,1) {mustBeReal} = 0
                w_val (1,1) {mustBeReal} = 0
            end
            
            obj.uf = u_val * ones(obj.nx, obj.ncy, obj.ncz);
            obj.vf = v_val * ones(obj.ncx, obj.ny, obj.ncz);
            if obj.dim == 3
                obj.wf = w_val * ones(obj.ncx, obj.ncy, obj.nz);
            else
                obj.wf = zeros(obj.ncx, obj.ncy, obj.nz);
            end
            
            fprintf('均匀速度场设置: u=%.3f, v=%.3f, w=%.3f\n', u_val, v_val, w_val);
        end


        function obj = set_scheme(obj, scheme)
            % 设置对流项离散格式
            % scheme: 'upwind' 或 'centraldifference'
            arguments
                obj
                scheme {mustBeMember(scheme, {'upwind', 'centraldifference'})}
            end
            obj.scheme = scheme;
            fprintf('对流项离散格式设置为: %s\n', scheme);
        end


        function obj = create_coefficient_data(obj)
            % 创建系数存储
            
            % 设置系数个数
            if obj.dim == 2
                obj.ncoef = 6; % aP, aE, aW, aN, aS, bsrc
            elseif obj.dim == 3
                obj.ncoef = 8; % aP, aE, aW, aN, aS, aT, aB, bsrc
            end
            
            % 系数存储位置索引
            obj.id_aP = 1;
            obj.id_aE = 2;
            obj.id_aW = 3;
            obj.id_aN = 4;
            obj.id_aS = 5;
            if obj.dim == 3
                obj.id_aT = 6;
                obj.id_aB = 7;
            end
            obj.id_bsrc = obj.ncoef;
            
            % 初始化系数数组
            obj.coef = zeros(obj.ncx, obj.ncy, obj.ncz, obj.ncoef);
        end
    
        % 设置非线性迭代次数
        function obj = set_nsteps(obj, nsteps)
            obj.nsteps = nsteps;  
        end
    
        % 设置时间步长
        function obj = set_dt(obj, dt)
            obj.dt = dt;
        end

        function obj = set_temperature_solver_param(obj, niter_t, relax_t, res_t, temp_tol)
            % 设置温度求解器参数
            % niter_t   - 温度线性求解器迭代次数
            % relax_t   - 松弛因子
            % res_t     - 线性求解器收敛误差
            % temp_tol  - 非线性迭代收敛误差
            
            arguments
                obj
                niter_t (1,1) {mustBePositive, mustBeInteger}
                relax_t (1,1) {mustBePositive, mustBeLessThan(relax_t, 2)}
                res_t (1,1) {mustBePositive}
                temp_tol (1,1) {mustBePositive}
            end
            
            obj.niter_t  = niter_t;
            obj.relax_t  = relax_t;
            obj.res_t    = res_t;
            obj.temp_tol = temp_tol;
        end
    end
end



