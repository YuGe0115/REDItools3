# 导入必要的模块
from reditools.plotter import Plotter

import matplotlib.pyplot as plt
import numpy as np

# 定义 ManhattanPlotter 类，继承自 Plotter 抽象基类
class ManhattanPlotter(Plotter):
    def plot(self):
        """
        实现曼哈顿图的绘图逻辑，覆盖 Plotter 类的抽象方法 plot。

        功能:
            - 按 'Region'（染色体）分组，绘制每个染色体的编辑事件散点图。
            - X 轴表示基因组位置（Position + 染色体偏移量），Y 轴表示编辑频率（Frequency）。
            - 不同染色体使用不同颜色，添加图例和标签。
            - 调用父类的 save 方法保存图表。

        输入:
            - self.data (pd.DataFrame): 包含 'Region', 'Position', 'Frequency' 字段的 Pandas DataFrame。
            - self.config (dict): 配置字典，包含绘图参数（如 figsize、colors、title）。
            - self.output_path (str): 图表保存路径（不含扩展名）。

        输出:
            - 生成曼哈顿图并保存到 self.output_path.<format>（默认 PNG 格式）。

        示例:
            - 如果 self.data 包含 'Region'（如 chr1, chr2）、'Position' 和 'Frequency' 列，
              self.config = {"figsize": (12, 6), "colors": ["#1f77b4", "#ff7f0e"]},
              self.output_path = "output/manhattan"，则生成曼哈顿图并保存为 output/manhattan.png。
        """
        plt.figure(figsize=self.config.get("figsize", (12, 6)))
        
        # 获取所有唯一的染色体（Region），如 ['chr1', 'chr2', 'chrM']
        chromosomes = self.data["Region"].unique()
        
        # 获取颜色列表，默认为蓝色、橙色、绿色、红色，循环使用
        # self.config.get("colors", [...]): 从配置获取颜色列表，默认值是 Matplotlib 的常用颜色
        colors = self.config.get("colors", ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"])
        
        # 初始化 X 轴偏移量，用于在 X 轴上分隔不同染色体
        x_offset = 0
        
        # 初始化列表，用于存储所有点的 X 轴位置（用于绘图）
        x_positions = []
        
        # 初始化列表，用于存储 X 轴刻度标签（染色体名称和位置）
        x_labels = []
        
        # 遍历每个染色体，绘制散点图
        for i, chrom in enumerate(chromosomes):
            # 筛选当前染色体的数据
            chrom_data = self.data[self.data["Region"] == chrom]
            
            # 计算当前染色体的 X 轴位置：基因组位置 + 偏移量
            # x_offset 确保不同染色体在 X 轴上不重叠
            x_pos = chrom_data["Position"] + x_offset
            
            # 将当前染色体的 X 轴位置添加到总列表
            x_positions.extend(x_pos)
            
            # 绘制散点图
            # x_pos: X 轴位置（Position + 偏移量）
            # chrom_data["Frequency"]: Y 轴值（编辑频率）
            # color=colors[i % len(colors)]: 循环使用颜色列表中的颜色
            # label=chrom: 为当前染色体设置图例标签
            # alpha=0.6: 设置点的透明度，防止重叠点过于密集
            # s=20: 设置点的大小
            plt.scatter(x_pos, chrom_data["Frequency"], color=colors[i % len(colors)], 
                        label=chrom, alpha=0.6, s=20)
            
            # 计算当前染色体的 X 轴刻度位置（使用 Position 的均值）和标签（染色体名称）
            x_labels.append((x_offset + chrom_data["Position"].mean(), chrom))
            
            # 更新偏移量：当前染色体的最大 Position + 1000000（染色体间距）
            # 1000000 是一个固定的间距值，确保染色体在 X 轴上分隔开
            x_offset += chrom_data["Position"].max() + 1000000
        
        # 设置图表标题
        # self.config.get("title", "RNA 编辑曼哈顿图"): 从 config 获取标题，默认为中文
        plt.title(self.config.get("title", "RNA Manhattan"))
        
        # 设置 X 轴标签为“基因组位置”
        plt.xlabel("Genomic position")
        
        # 设置 Y 轴标签为“编辑频率”
        plt.ylabel("Frequency")
        
        # 设置 X 轴刻度
        # [x for x, _ in x_labels]: 提取 x_labels 中的位置（均值 Position + 偏移量）
        # [chrom for _, chrom in x_labels]: 提取 x_labels 中的染色体名称
        # rotation=45: 旋转 X 轴标签 45 度，便于阅读
        plt.xticks([x for x, _ in x_labels], [chrom for _, chrom in x_labels], rotation=45)
        
        # 添加图例，显示每个染色体的颜色和名称
        plt.legend()
        
        # 调用父类的 save 方法保存图表
        # self.config.get("format", "png"): 从 config 获取保存格式，默认 "png"
        self.save(self.config.get("format", "png"))