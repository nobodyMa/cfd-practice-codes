# cfd-practice-codes
这些是我在学习CFD过程中，编写的玩具代码
## 项目名称
### 一、模仿的代码
#### 1.heat_conduction
##### 参考
| 类型 | 名称 | 来源/链接 | 说明 |
|------|------|------------|------|
| 代码参考 | cfd0to1 | [CFD从0到1 Lecture 1.40 代码分享总述，代码](https://www.bilibili.com/video/BV1Zs4y1Q7dk/?spm_id_from=333.1387.collection.video_card.click&vd_source=0745441b4a83ceba73d32af3b7b0a955) | 参考程序的整体逻辑 |
| 书籍 | The Finite Volume Method in Computational Fluid Dynamics (F. Moukalled  L. Mangani  M. Darwish) |  | 中心差分法和一阶迎风格式 |
| 数据 |  | [2D conduction in rectangular plate](https://www.theseus-fe.com/resources/validations/basic-heat-transfer) |x=0.5和x=0.0833两个位置的参考温度分布|
##### 开发时间
- **2025.10.30 - 2025.11.02**：完成程序编写，但是该版本没有上传Github，已经丢失。
- **22025.11.08**：修改雅可比迭代推出的判断逻辑，并添加一阶迎风格式。

#### 2.heat_convection
##### 参考
| 类型 | 名称 | 来源/链接 | 说明 |
|------|------|------------|------|
| 代码参考 | cfd0to1 | [CFD从0到1 Lecture 1.40 代码分享总述，代码](https://www.bilibili.com/video/BV1Zs4y1Q7dk/?spm_id_from=333.1387.collection.video_card.click&vd_source=0745441b4a83ceba73d32af3b7b0a955) | 参考程序的整体逻辑 |
| 书籍 | The Finite Volume Method in Computational Fluid Dynamics (F. Moukalled  L. Mangani  M. Darwish) |  | 热对流的离散方程和一维对流扩散问题的解析解 |
##### 开发时间
- **2025.11.09**：在heat_conduction的基础上添加热对流的求解，并复现The Finite Volume Method in Computational Fluid Dynamics (F. Moukalled  L. Mangani  M. Darwish)上的fig11.7一维对流扩散问题中心点无量纲温度随Pe数变化。


### 二、复现的代码

