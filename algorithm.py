import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.rcParams['font.sans-serif']=['SimHei']#显示中文
plt.rcParams['axes.unicode_minus']=False#显示负号

def desmoother(surface,num=100):
    m,n = surface.shape
    for _ in range(num):
        size  = np.random.choice(np.arange(10,31,1))
        scale = np.random.choice([0.2,0.3,0.5,0.7,0.8])
        k = (size+1)//2
        patch = np.ones(shape=(size,size))
        for i in range(k):
            patch[i+1:size-(i+1),i+1:size-(i+1)] = 1+0.1*(i+1)
        start_x = np.random.choice(np.arange(0,n-size,1))
        start_y = np.random.choice(np.arange(0,m-size,1))

        surface[start_y:start_y+size,start_x:start_x+size] += scale*patch
        
    return surface

def get_point():
    #平面、斜坡、凹地、凸地
    names = ['平面','斜坡','凹地','凸地']
    X = np.array([-10,10,-10,10,0])
    Y = np.array([10,10,-10,-10,0])
    Z = np.array([[0,0,0,0,0],[-5,5,-5,5,0],[2,2,2,2,0],[-2,-2,-2,-2,0]])
    for i in range(4):
        yield X,Y,Z[i],names[i]
    
def surface_fitting(X,Y,Z):

    def f(x): #处理打印方程时的符号问题

        round(x,2) #保留两位小数

        return '+'+str(x) if x>= 0 else str(x)

    n = len(X)
    
    #求方程系数
    sigma_x1 = np.sum(X)
    sigma_y1 = np.sum(Y)
    sigma_z1 = np.sum(Z)

    sigma_x2 = np.sum(X**2)
    sigma_y2 = np.sum(Y**2)

    sigma_x3 = np.sum(X**3)
    sigma_y3 = np.sum(Y**3)

    sigma_x4 = np.sum(X**4)
    sigma_y4 = np.sum(Y**4)

    sigma_x1y1 = np.sum(X*Y)
    sigma_x1y2 = np.sum(X*Y**2)
    sigma_x1y3 = np.sum(X*Y**3)
    sigma_x2y1 = np.sum(X**2*Y)
    sigma_x3y1 = np.sum(X**3*Y)
    sigma_x2y2 = np.sum(X**2*Y**2)

    sigma_z1x1 = np.sum(Z*X)
    sigma_z1y1 = np.sum(Z*Y)
    sigma_z1x2 = np.sum(Z*X**2)
    sigma_z1y2 = np.sum(Z*Y**2)

    sigma_z1x1y1 = np.sum(X*Y*Z)

    #给出对应方程的矩阵形式                        
    a = np.array([[sigma_x4,sigma_x2y2,sigma_x3y1,sigma_x3,sigma_x2y1,sigma_x2],

                 [sigma_x2y2,sigma_y4,sigma_x1y3,sigma_x1y2,sigma_y4,sigma_y2],

                 [sigma_x3y1,sigma_x1y3,sigma_x2y2,sigma_x2y1,sigma_x1y2,sigma_x1y1],

                 [sigma_x3,sigma_x1y2,sigma_x2y1,sigma_x2,sigma_x1y1,sigma_x1],

                 [sigma_x2y1,sigma_y4,sigma_x1y2,sigma_x1y1,sigma_y2,sigma_y1],

                 [sigma_x2,sigma_y2,sigma_x1y1,sigma_x1,sigma_y1,n]])

    b = np.array([sigma_z1x2,sigma_z1y2,sigma_z1x1y1,sigma_z1x1,sigma_z1y1,sigma_z1])

    #判断是否满足唯一解条件，如果满足则用高斯消元解线性方程
    #print(np.linalg.matrix_rank(a))
    #print(np.linalg.matrix_rank(np.concatenate((a,b.reshape(6,1)),axis=1)))
    if np.linalg.matrix_rank(a) == np.linalg.matrix_rank(np.concatenate((a,b.reshape(6,1)),axis=1)) == 6:
        res = np.linalg.solve(a,b)     #ax = b , res <- x
    else:
        raise Exception('输入的观测点无效，请重新输入！')
    
    #输出方程形式(%.6s表示保留6个字符)
    print("z = %.6s*x^2%.6s*y^2%.6s*xy%.6s*x%.6s*y%.6s"%(f(res[0]),f((res[1]+0.9*res[0])),f(res[2]),f(res[3]),f(res[4]),f(res[5])))#(res[1]+0.99*res[0])是为了修复算法拟合凹面或凸面时有一侧不变的BUG

    return res

#主函数
if __name__ == '__main__':
    
    fig = plt.figure(dpi=128)#建立一个空间
    ax = Axes3D(fig)# 3D坐标
    u = np.linspace(-20,20,256) # 创建一个等差数列
    x, y = np.meshgrid(u,u) #转化成矩阵
    point_generator = get_point()
    
    for _ in range(4):
        
        #清除原有图像
        plt.cla()
        
        #仪器每次输入五个观测点
        #X,Y,Z,name = next(point_generator)
        
        X = np.array([-10, 10, -10, 10, 0])
        Y = np.array([10, 10, -10, -10, 0])
        Z = np.array([1, 1, 1, 1, 0])
        name = '凹面'
        
        
        #得到拟合曲面方程
        print(name,":")
        res = surface_fitting(X,Y,Z)
        z = res[0]*x**2+(res[1]+0.99*res[0])*y**2+res[2]*x*y+res[3]*x+res[4]*y+res[5]#(res[1]+0.99*res[0])是为了修复算法拟合凹面或凸面时有一侧不变的BUG
        #加入噪声使曲面更真实
        z += 0.6*np.exp(np.sin(0.6*np.cos(x+y)))-0.3*np.exp(np.sin(0.3*np.exp(np.sin(x))-3*np.sin(y))) + 0.1*z*np.sin(0.01*x*y) + 0.6*res[0]*y**2
        #随机加入凹凸块
        z = desmoother(z,num=100)
        
        #实时绘制曲面图和离散点
        ax.plot_surface(x,y,z,rstride=3,cstride=3)#cmap=plt.get_cmap('cool')
        #ax.scatter(X, Y, Z,c='black',s=100)
        ax.xaxis.axes.set_xlim3d(left=-20,right=20)#设置x轴范围
        ax.yaxis.axes.set_ylim3d(bottom=-20,top=20)#设置y轴范围
        ax.zaxis.axes.set_zlim3d(bottom=-20,top=20)#设置z轴范围
        
        plt.axis('off')
        plt.title(name)
        plt.pause(1)
        #break
    
