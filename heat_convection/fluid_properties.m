classdef fluid_properties
    properties
        % 属性场
        mesh
        dens      % 密度场
        mu        % 粘度场 
        con       % 传导系数
        spheat    % 比热

        % 热源
        heat_source
    end

    methods
        function obj = fluid_properties(mesh)
            % 初始化流体属性
            % mesh_obj - structured_mesh 对象
            arguments
                mesh (1,1) structured_mesh
            end
            obj.mesh = mesh;
            % 从网格对象获取尺寸
            ncx = obj.mesh.ncx;
            ncy = obj.mesh.ncy;
            ncz = obj.mesh.ncz;
            
            % 初始化场数据
            obj.dens = ones(ncx, ncy, ncz);      % 密度场
            obj.mu   = ones(ncx, ncy, ncz);      % 粘度场 

            % 材料属性
            obj.con   = 0.0;            % 传导系数
            obj.spheat  = 0.0;            % 比热

            % 热源
            obj.heat_source = 0.0;
        end
        
        function obj = set_initial_fluid(obj, dens, mu, con_val, spheat_val, heat_source_val)
            % 自定义密度和粘度场
            % dens_val - 密度值
            % mu_val   - 粘度值
            % con_val  - 传导系数
            % spheat_val - 比热
            % heat_source_val - 热源

            arguments
                obj
                dens (1,1) {mustBeReal, mustBePositive}
                mu (1,1) {mustBeReal, mustBePositive}
                con_val (1,1) {mustBeReal, mustBeNonnegative}
                spheat_val (1,1) {mustBeReal, mustBeNonnegative}
                heat_source_val  (1,1) {mustBeReal} = 0
            end
            
            % 从网格对象获取尺寸
            ncx = obj.mesh.ncx;
            ncy = obj.mesh.ncy;
            ncz = obj.mesh.ncz;
            
            % 设置密度和粘度场
            obj.dens = dens * ones(ncx, ncy, ncz);    % 密度场
            obj.mu   = mu * ones(ncx, ncy, ncz);    % 粘度场
            obj.con   = con_val;        % 传导系数
            obj.spheat  = spheat_val;       % 比热
            obj.heat_source = heat_source_val;

            % display_info(obj);
        end
        
        function display_info(obj)
            % 显示流体信息
            fprintf('=== 流体属性信息 ===\n');
            fprintf('网格尺寸: %d × %d × %d\n', obj.mesh.ncx, obj.mesh.ncy, obj.mesh.ncz);
            fprintf('密度场范围: [%.3f, %.3f]\n', min(obj.dens(:)), max(obj.dens(:)));
            fprintf('粘度场范围: [%.3f, %.3f]\n', min(obj.mu(:)), max(obj.mu(:)));
            fprintf('传导系数: %.6f\n', obj.con);
            fprintf('比热: %.6f\n', obj.spheat);
            fprintf('热源: %.6f\n', obj.heat_source);
        end

    end
end