## 1 introduction
多智能体编队引起广泛关注，同时反无人机技术的社会需求也日益增加。但是多智能体编队在反无人机中的应用相关的工作还比较少。
阐述编队控制方法，论证角度编队的优越性。

## 2 Problem formulation

使用双积分模型 绳索张力模型

#### Assumption1: 
所有的agent配备LOS，可以测量
可以获得编队控制所需要的相对bearing信息。

#### Assumption2: 
Agent1可以测量相邻的agent的距离信息，并且可以通过角度信息计算对边距离大小，并且将距离信息通信给相邻agent。
可以控制编队缩放；并且使得agent1 2 4 获取在其局部坐标系下的平面中心位置，进而控制平面法向量。

#### Assumption3:
初始位置在期望构型附近
可以

引出角度定义

$$\measuredangle ijk=\arccos\left(b_{ij}^\top b_{ik}\right)$$

其中$\measuredangle jik\in[0,\pi]$ 且 $\measuredangle ijk = \measuredangle kij$

引出共面项定义

引出法向量误差定义(包括分支)

关于无向角度的说明，即表示无向角度在当前的捕网场景需求不符合 remark?

defination1: 3D空间中的角度刚性 是否引用

defination2: 通过图论

阐述角度刚性

#### Problem1

......->0
......->0
......->0
......->0
......->0

#### Theorem1
当前四个agent情况下当一个agent离开指定参考平面，角误差和共面误差 proof:(一阶误差和二阶误差的方式说明原来的角度约束对于离面误差不敏感):

#### remark1
解释为什么要加共面项

## 3 Main Resluts

#### Theorem1
当前约束存在冗余 proof:

#### Theorem2 
当前控制律在邻域内收敛 proof:


