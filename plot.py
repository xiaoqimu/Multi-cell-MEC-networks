import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import  ConnectionPatch
import warnings
warnings.filterwarnings('ignore')
from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import FuncFormatter
import random
data0 = np.loadtxt(open("result0_0.csv","rb"),delimiter=",",skiprows=0) 
plt.figure()
plt.figure(figsize=(7,7))#将画布设定为正方形，则绘制的饼图是正圆
label=['Local','Beijing','Chengdu','Singapore','Sydney','Silicon Valley']#定义饼图的标签，标签是列表
explode=[0.0,0.0,0.0,0.8,0.05,0.45]#设定各项距离圆心n个半径
values=[data0[0],data0[1],data0[3],data0[4],data0[5],data0[6]]
colors = ['#8ECFC9', '#FFBE7A', '#FA7F6F', '#82B0D2', '#BEB8DC','#E7DAD2']

plt.pie(values,explode=explode,colors=colors, autopct='%1.4f%%',pctdistance=0.6,textprops={'fontsize': 12})#绘制饼图
plt.legend(labels=label,loc='upper left')
plt.savefig('preference')#保存图片
########################在线################
data = np.loadtxt(open("result9_21.csv","rb"),delimiter=",",skiprows=0) 
plt.figure()
l=len(data[0][:])
lh=l//50
x=np.arange(0, l, lh)
x1=np.arange(0, l, 5)
line1,=plt.plot(x,data[0][:][0:l:lh],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line2,=plt.plot(x,data[5][:][0:l:lh],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line3,=plt.plot(x,data[10][:][0:l:lh],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line4,=plt.plot(x,data[15][:][0:l:lh],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line5,=plt.plot(x,data[20][:][0:l:lh],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
ll=plt.legend([line1,line2,line3,line4,line5],["DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","D-DEBO"],loc='upper left')#添加图例
plt.ylabel("Accumulated regret",fontsize=14)#设置纵轴单位
plt.xlabel("Time slot t",fontsize=14)#设置横轴单位
plt.savefig('Accumulated_regret.png')
plt.figure()

line1,=plt.plot(x,data[1][:][0:l:lh],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line2,=plt.plot(x,data[6][:][0:l:lh],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line3,=plt.plot(x,data[11][:][0:l:lh],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line4,=plt.plot(x,data[16][:][0:l:lh],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line5,=plt.plot(x,data[21][:][0:l:lh],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
ll=plt.legend([line1,line2,line3,line4,line5],["DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","D-DEBO"],loc='upper left')#添加图例
plt.ylabel("Time average regret",fontsize=14)#设置纵轴单位
plt.xlabel("Time slot t",fontsize=14)#设置横轴单位
plt.savefig('Time_average_regret.png')

plt.figure()

line1,=plt.plot(x,data[2][:][0:l:lh],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line2,=plt.plot(x,data[7][:][0:l:lh],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line3,=plt.plot(x,data[12][:][0:l:lh],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line4,=plt.plot(x,data[17][:][0:l:lh],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line5,=plt.plot(x,data[22][:][0:l:lh],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
ll=plt.legend([line1,line2,line3,line4,line5],["DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","D-DEBO"],loc='upper left')#添加图例
plt.ylabel("Accumulated reward",fontsize=14)#设置纵轴单位
plt.xlabel("Time slot t",fontsize=14)#设置横轴单位
plt.savefig('Accumulated_reward.png')

plt.figure()

line1,=plt.plot(x,data[3][:][0:l:lh],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line2,=plt.plot(x,data[8][:][0:l:lh],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line3,=plt.plot(x,data[13][:][0:l:lh],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line4,=plt.plot(x,data[18][:][0:l:lh],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line5,=plt.plot(x,data[23][:][0:l:lh],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
ll=plt.legend([line1,line2,line3,line4,line5],["DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","D-DEBO"],loc='upper left')#添加图例
plt.ylabel("Time average reward",fontsize=14)#设置纵轴单位
plt.xlabel("Time slot t",fontsize=14)#设置横轴单位
plt.savefig('Time_average_reward.png')

plt.figure()

line1,=plt.plot(x,data[4][:][0:l:lh],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line2,=plt.plot(x,data[9][:][0:l:lh],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line3,=plt.plot(x,data[14][:][0:l:lh],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line4,=plt.plot(x,data[19][:][0:l:lh],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line5,=plt.plot(x,data[24][:][0:l:lh],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
ll=plt.legend([line1,line2,line3,line4,line5],["DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","D-DEBO"],loc='upper left')#添加图例
plt.ylabel("Reward ratio",fontsize=14)#设置纵轴单位
plt.xlabel("Time slot t",fontsize=14)#设置横轴单位
plt.savefig('Reward_ratio.png')

####################################################
def zone_and_linked(ax,axins,zone_left,zone_right,x,y,index_ls,linked,
                    x_ratio,y_ratio):
    """缩放内嵌图形，并且进行连线
    ax:         调用plt.subplots返回的画布。例如 fig,ax = plt.subplots(1,1)
    axins:      内嵌图的画布。 例如 axins = ax.inset_axes((0.4,0.1,0.4,0.3))
    zone_left:  要放大区域的横坐标左端点
    zone_right: 要放大区域的横坐标右端点
    x:          X轴标签
    y:          列表所有y值
    linked:     进行连线的位置，{'bottom','top','left','right'}
    x_ratio:    X轴缩放比例
    y_ratio:    Y轴缩放比例
    """
    xlim_left = x[zone_left]-(x[zone_right]-x[zone_left])*x_ratio
    xlim_right = x[zone_right]+(x[zone_right]-x[zone_left])*x_ratio

    y_data = np.hstack([yi[zone_left:zone_right] for yi in y])
    ylim_bottom = np.min(y_data)-(np.max(y_data)-np.min(y_data))*y_ratio
    ylim_top = np.max(y_data)+(np.max(y_data)-np.min(y_data))*y_ratio

    axins.set_xlim(xlim_left-0.1, xlim_right+0.1)
    axins.set_ylim(ylim_bottom, ylim_top)

    ax.plot([xlim_left,xlim_right,xlim_right,xlim_left,xlim_left],[ylim_bottom,ylim_bottom,ylim_top,ylim_top,ylim_bottom],"black")
    #plt.xticks([])  # 去掉横坐标值
    if linked=='right':
        xyA_1, xyB_1 = (xlim_left,ylim_bottom), (xlim_left,ylim_top)
        xyA_2, xyB_2 = (xlim_right,ylim_bottom), (xlim_right,ylim_top)
    else:
        xyA_1, xyB_1 = (xlim_left,ylim_top), (xlim_left,ylim_bottom)
        xyA_2, xyB_2 = (xlim_right,ylim_top), (xlim_right,ylim_bottom)
    con = ConnectionPatch(xyA=xyA_1,xyB=xyB_1,coordsA="data",
                          coordsB="data",axesA=axins,axesB=ax)
    axins.add_artist(con)
    con = ConnectionPatch(xyA=xyA_2,xyB=xyB_2,coordsA="data",
                          coordsB="data",axesA=axins,axesB=ax)
    axins.add_artist(con)

def zone_and_linked1(ax,axins,zone_left,zone_right,x,y,index_ls,linked,
                    x_ratio,y_ratio):
    """缩放内嵌图形，并且进行连线
    ax:         调用plt.subplots返回的画布。例如 fig,ax = plt.subplots(1,1)
    axins:      内嵌图的画布。 例如 axins = ax.inset_axes((0.4,0.1,0.4,0.3))
    zone_left:  要放大区域的横坐标左端点
    zone_right: 要放大区域的横坐标右端点
    x:          X轴标签
    y:          列表所有y值
    linked:     进行连线的位置，{'bottom','top','left','right'}
    x_ratio:    X轴缩放比例
    y_ratio:    Y轴缩放比例
    """
    xlim_left = x[zone_left]-(x[zone_right]-x[zone_left])*x_ratio
    xlim_right = x[zone_right]+(x[zone_right]-x[zone_left])*x_ratio

    y_data = np.hstack([yi[zone_left:zone_right] for yi in y])
    ylim_bottom = np.min(y_data)-(np.max(y_data)-np.min(y_data))*y_ratio
    ylim_top = np.max(y_data)+(np.max(y_data)-np.min(y_data))*y_ratio

    axins.set_xlim(xlim_left-0.1, xlim_right-0.9)
    axins.set_ylim(ylim_bottom, ylim_top)

    #ax.plot([xlim_left,xlim_right,xlim_right,xlim_left,xlim_left],[ylim_bottom,ylim_bottom,ylim_top,ylim_top,ylim_bottom],"black")
    #plt.xticks([])  # 去掉横坐标值
    if linked=='right':
        xyA_1, xyB_1 = (xlim_left,ylim_bottom), (xlim_left,ylim_top)
        xyA_2, xyB_2 = (xlim_right,ylim_bottom), (xlim_right,ylim_top)
    else:
        xyA_1, xyB_1 = (xlim_left,ylim_top), (xlim_left,ylim_bottom)
        xyA_2, xyB_2 = (xlim_right,ylim_top), (xlim_right,ylim_bottom)
    con = ConnectionPatch(xyA=xyA_1,xyB=xyB_1,coordsA="data",
                          coordsB="data",axesA=axins,axesB=ax)
    axins.add_artist(con)
   # con = ConnectionPatch(xyA=xyA_2,xyB=xyB_2,coordsA="data",
                          #coordsB="data",axesA=axins,axesB=ax)
   # axins.add_artist(con)


data1 = np.loadtxt(open("result1_89_u.csv","rb"),delimiter=",",skiprows=0) 
index_ls=['i','ii','iii','iv','v','vi']
scale_ls = range(6)
plt.figure()
fig, ax = plt.subplots(1,1)
line1,=ax.plot(data1[0][:],color = '#c82423', marker='*', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)
line2,=ax.plot(data1[1][:],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line3,=ax.plot(data1[2][:],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line4,=ax.plot(data1[3][:],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line5,=ax.plot(data1[4][:],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line6,=ax.plot(data1[5][:],color = '#32B897',marker='X', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line7,=ax.plot(data1[6][:],color = '#999999',marker='d', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line8,=ax.plot(data1[7][:],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上

ax.legend([line1,line2,line3,line4,line5,line6,line7,line8],["DGOSS","DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","J-UACO","PSO","D-DEBO"],loc='upper left')
plt.ylabel("Time average QoE",fontsize=14)#设置纵轴单位
plt.xlabel("Network size",fontsize=14)#设置横轴单位
# 绘制缩放图

axins = ax.inset_axes((0.7, 0.09, 0.2, 0.15))
# 在缩放图中也绘制主图所有内容，然后根据限制横纵坐标来达成局部显示的目的
axins.plot(scale_ls,data1[1][:],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
axins.plot(scale_ls,data1[2][:],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[3][:],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[4][:],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[7][:],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上


# 局部显示并且进行连线
zone_and_linked1(ax, axins, 0, 1, scale_ls , [data1[1][:],data1[7][:]],index_ls, 'right',0,0.3)
plt.xticks(scale_ls,index_ls,fontsize=14) 
plt.savefig('QoE.png')

plt.figure()
# 改变文字大小参数-fontsize
labels =['i','ii','iii','iv','v','vi']
x = np.arange(len(labels))  # 标签位置
width = 0.2  # 柱状图的宽度"DGOSS","DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","J-UACO","PSO"
rects1 = plt.bar(x -width, data1[0][:], width, label='DGOSS', hatch="...", color='#c82423', edgecolor="k")
rects2 =  plt.bar(x, data1[5][:], width, label='J-UACO', hatch="***", color='#32B897', edgecolor="k")
rects3 =  plt.bar(x +width, data1[6][:], width, label='PSO', hatch="xxx", color='#999999', edgecolor="k")
# 为y轴、标题和x轴等添加一些文本。
plt.ylabel("Time average QoE",fontsize=14)#设置纵轴单位
plt.xlabel("Network size",fontsize=14)#设置横轴单位
plt.xticks(x,labels,fontsize=14) 
plt.legend(loc='upper left')
plt.savefig('QoE_offline.png')

plt.figure()
# 绘制主图
scale_ls = [1,2,3,4,5,6]
fig, ax = plt.subplots(1,1)
line1,=ax.plot(scale_ls,data1[8][:],color = '#c82423', marker='*', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)
line2,=ax.plot(scale_ls,data1[9][:],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line3,=ax.plot(scale_ls,data1[10][:],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line4,=ax.plot(scale_ls,data1[11][:],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line5,=ax.plot(scale_ls,data1[12][:],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line6,=ax.plot(scale_ls,data1[13][:],color = '#32B897',marker='X', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line7,=ax.plot(scale_ls,data1[14][:],color = '#999999',marker='d', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line8,=ax.plot(scale_ls,data1[15][:],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上

ax.legend([line1,line2,line3,line4,line5,line6,line7,line8],["DGOSS","DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","J-UACO","PSO","D-DEBO"],loc='upper left')
plt.ylabel("Time average running time (ms)",fontsize=14)#设置纵轴单位
plt.xlabel("Network size",fontsize=14)#设置横轴单位
# 绘制缩放图
axins = ax.inset_axes((0.4, 0.7, 0.3, 0.2))
# 在缩放图中也绘制主图所有内容，然后根据限制横纵坐标来达成局部显示的目的
axins.plot(scale_ls,data1[9][:],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
axins.plot(scale_ls,data1[10][:],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[11][:],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[12][:],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[15][:],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
# 局部显示并且进行连线
zone_and_linked(ax, axins, 3, 5, scale_ls , [data1[9][:],data1[10][:],data1[11][:],data1[12][:],data1[15][:]],index_ls, 'right',0,0.2)

plt.xticks(scale_ls,index_ls,fontsize=14) 

plt.savefig('Running_time.png')


plt.figure()
# 绘制主图
scale_ls = [1,2,3,4,5,6]
fig, ax = plt.subplots(1,1)
line2,=ax.plot(scale_ls,data1[9][:],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line3,=ax.plot(scale_ls,data1[10][:],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line4,=ax.plot(scale_ls,data1[11][:],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line5,=ax.plot(scale_ls,data1[12][:],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line8,=ax.plot(scale_ls,data1[15][:],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上

ax.legend([line1,line2,line3,line4,line5,line6,line7,line8],["DGOSS","DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","J-UACO","PSO","D-DEBO"],loc='upper left')
plt.ylabel("Time average running time (ms)",fontsize=14)#设置纵轴单位
plt.xlabel("Network size",fontsize=14)#设置横轴单位
# 绘制缩放图
plt.xticks(scale_ls,index_ls) 
plt.xticks(scale_ls,index_ls) 
plt.savefig('Running_time_online.png')


data1 = np.loadtxt(open("result9_11.csv","rb"),delimiter=",",skiprows=0) 
index_ls=['1','2','3','4','5','6','7','8','9','10']
scale_ls = range(10)
plt.figure()
line1,=plt.plot(data1[0][:],color = '#c82423', marker='*', markerfacecolor='none',markeredgewidth='1',linewidth = '1')
line2,=plt.plot(data1[1][:],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line3,=plt.plot(data1[2][:],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line4,=plt.plot(data1[3][:],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line5,=plt.plot(data1[4][:],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line6,=plt.plot(data1[5][:],color = '#32B897',marker='X', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line7,=plt.plot(data1[6][:],color = '#999999',marker='d', markerfacecolor='none',markeredgewidth='1',linewidth = '1')#同上
line8,=plt.plot(data1[7][:],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
plt.xticks(scale_ls,index_ls,fontsize=14) 
ll=plt.legend([line1,line2,line3,line4,line5,line6,line7,line8],["DGOSS","DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","J-UACO","PSO","D-DEBO"],loc='upper left')#添加图例
plt.ylabel("Time average QoE",fontsize=14)#设置纵轴单位
plt.xlabel("Epoch",fontsize=14)#设置横轴单位
plt.savefig('QoE_2.png')

plt.figure()
# 改变文字大小参数-fontsize
x = np.arange(len(index_ls))  # 标签位置
width = 0.2  # 柱状图的宽度"DGOSS","DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","J-UACO","PSO"
rects1 = plt.bar(x -width, data1[0][:], width, label='DGOSS', hatch="...", color='#c82423', edgecolor="k")
rects2 =  plt.bar(x, data1[5][:], width, label='J-UACO', hatch="***", color='#32B897', edgecolor="k")
rects3 =  plt.bar(x +width, data1[6][:], width, label='PSO', hatch="xxx", color='#999999', edgecolor="k")
# 为y轴、标题和x轴等添加一些文本。
plt.ylabel("Time average QoE",fontsize=14)#设置纵轴单位
plt.xlabel("Epoch",fontsize=14)#设置横轴单位
plt.xticks(x,index_ls,fontsize=14) 
plt.legend(loc='upper left')
plt.savefig('QoE2_offline.png')


plt.figure()
# 绘制主图
scale_ls = [1,2,3,4,5,6,7,8,9,10]
fig, ax = plt.subplots(1,1)
line1,=ax.plot(scale_ls,data1[8][:],color = '#c82423', marker='*', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)
line2,=ax.plot(scale_ls,data1[9][:],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
line3,=ax.plot(scale_ls,data1[10][:],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line4,=ax.plot(scale_ls,data1[11][:],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line5,=ax.plot(scale_ls,data1[12][:],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line6,=ax.plot(scale_ls,data1[13][:],color = '#32B897',marker='X', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line7,=ax.plot(scale_ls,data1[14][:],color = '#999999',marker='d', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
line8,=ax.plot(scale_ls,data1[15][:],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上

ax.legend([line1,line2,line3,line4,line5,line6,line7,line8],["DGOSS","DTE-DOL", "M-UCB","M-Exp3","M-Exp3-IX","J-UACO","PSO","D-DEBO"],loc='upper left')
plt.ylabel("Time average running time (ms)",fontsize=14)#设置纵轴单位
plt.xlabel("Epoch",fontsize=14)#设置横轴单位
# 绘制缩放图
axins = ax.inset_axes((0.3, 0.14, 0.6, 0.15))
# 在缩放图中也绘制主图所有内容，然后根据限制横纵坐标来达成局部显示的目的
axins.plot(scale_ls,data1[8][:],color = '#c82423', marker='*', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)
axins.plot(scale_ls,data1[9][:],color = '#FFBE7A', marker='o', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#给曲线设置标识。并把曲线赋给一个变量，方便下面添加图例时候应用
axins.plot(scale_ls,data1[10][:],color = '#BEB8DC',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[11][:],color = '#82B0D2',marker='s', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[12][:],color = '#FA7F6F',marker='D', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[13][:],color = '#32B897',marker='X', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[14][:],color = '#999999',marker='d', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上
axins.plot(scale_ls,data1[15][:],color = '#32B897',marker='^', markerfacecolor='none',markeredgewidth='1',linewidth = '1',alpha=0.7)#同上

# 局部显示并且进行连线
zone_and_linked(ax, axins, 4, 8, scale_ls , [data1[9][:],data1[10][:],data1[11][:],data1[12][:],data1[15][:]],index_ls, 'right',0,0.1)
plt.xticks(scale_ls,index_ls,fontsize=14) 
plt.savefig('zitu.png')

def to_percent(temp, position):
    return '%1.0f'%(100*temp) + '%'

# 分辨率参数-dpi，画布大小参数-figsize
plt.figure(figsize=(10,7))
# 改变文字大小参数-fontsize
labels = ['Beijing','Shanghai','Chengdu','Singapore','Sydney','Silicon Valley','London']
a =np.loadtxt(open("result0_1.csv","rb"),delimiter=",",skiprows=0)
a=np.transpose(a)
x = np.arange(len(labels))  # 标签位置
width = 0.1  # 柱状图的宽度
'/', '', '|', '-', '+', 'x', 'o', 'O', '.', '*'
rects1 = plt.bar(x - width * 4, a[0][:], width, label='Local', hatch="...", color='#c82423', edgecolor="k")
rects2 =  plt.bar(x - width*3, a[1][:], width, label='Beijing', hatch="ooo", color='#FFBE7A', edgecolor="k")
rects3 =  plt.bar(x -width*2, a[2][:], width, label='Shanghai', hatch="+++", color='#BEB8DC', edgecolor="k")
rects4 = plt.bar(x -width, a[3][:], width, label='Chengdu', hatch="xxx", color='#82B0D2', edgecolor="k")
rects5 =  plt.bar(x, a[4][:], width, label='Singapore', hatch="***", color='#32B897', edgecolor="k")
rects6 =  plt.bar(x + width ,a[5][:], width, label='Sydney', hatch="---", color='#FA7F6F', edgecolor="k")
rects7 = plt.bar(x + width*2, a[6][:], width, label='Silicon Valley', hatch="///", color='#999999', edgecolor="k")
rects8 = plt.bar(x + width * 3 , a[7][:], width, label='London', hatch="", color='w', edgecolor="k")
# 为y轴、标题和x轴等添加一些文本。
plt.ylabel("Task offloading decision ratio",fontsize=14)#设置纵轴单位
plt.xlabel("Cell",fontsize=14)#设置横轴单位
plt.xticks(x,labels,fontsize=14) 
plt.yticks(fontsize=14) 
plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
plt.legend(loc='upper left',fontsize=14)
plt.savefig('zhu.png')


plt.figure(figsize=(10,7))
# 改变文字大小参数-fontsize
labels = ['Beijing','Shanghai','Chengdu','Singapore','Sydney','Silicon Valley','London']
a =np.loadtxt(open("result1_1.csv","rb"),delimiter=",",skiprows=0)
a=np.transpose(a)
x = np.arange(len(labels))  # 标签位置
width = 0.1  # 柱状图的宽度
'/', '', '|', '-', '+', 'x', 'o', 'O', '.', '*'
rects1 = plt.bar(x - width * 4, a[0][:], width, label='Local', hatch="...", color='#c82423', edgecolor="k")
rects2 =  plt.bar(x - width*3, a[1][:], width, label='Beijing', hatch="ooo", color='#FFBE7A', edgecolor="k")
rects3 =  plt.bar(x -width*2, a[2][:], width, label='Shanghai', hatch="+++", color='#BEB8DC', edgecolor="k")
rects4 = plt.bar(x -width, a[3][:], width, label='Chengdu', hatch="xxx", color='#82B0D2', edgecolor="k")
rects5 =  plt.bar(x, a[4][:], width, label='Singapore', hatch="***", color='#32B897', edgecolor="k")
rects6 =  plt.bar(x + width ,a[5][:], width, label='Sydney', hatch="---", color='#FA7F6F', edgecolor="k")
rects7 = plt.bar(x + width*2, a[6][:], width, label='Silicon Valley', hatch="///", color='#999999', edgecolor="k")
rects8 = plt.bar(x + width * 3 , a[7][:], width, label='London', hatch="", color='w', edgecolor="k")
# 为y轴、标题和x轴等添加一些文本。
plt.ylabel("Task offloading decision ratio",fontsize=14)#设置纵轴单位
plt.xlabel("Cell",fontsize=14)#设置横轴单位
plt.xticks(x,labels,fontsize=14) 
plt.yticks(fontsize=14) 
plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
plt.legend(loc='upper left',fontsize=14)
plt.savefig('zhu1.png')



 

 
