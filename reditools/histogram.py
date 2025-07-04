# 导入必要的模块
from reditools.plotter import Plotter  # 从当前包的 plotter 模块导入 Plotter 抽象基类
import matplotlib.pyplot as plt
import seaborn as sns

class HistogramPlotter(Plotter):
    def plot(self):
        """
        实现直方图的绘图逻辑，覆盖 Plotter 类的抽象方法 plot。

        功能:
            - 使用 Seaborn 的 histplot 函数绘制直方图，展示 'Frequency' 字段的分布。
            - 从 self.config 获取用户自定义参数（如图像大小、标题、bin 数量）。
            - 设置图表的标题、X 轴和 Y 轴标签。
            - 调用父类的 save 方法保存图表。

        输入:
            - self.data (pd.DataFrame): 包含 'Frequency' 字段的 Pandas DataFrame。
            - self.config (dict): 配置字典，包含绘图参数（如 figsize、title、bins）。
            - self.output_path (str): 图表保存路径（不含扩展名）。

        输出:
            - 生成直方图并保存到 self.output_path.<format>（默认 PNG 格式）。

        示例:
            - 如果 self.data 包含 'Frequency' 列，self.config = {"figsize": (12, 6), "bins": 50}，
              self.output_path = "output/histogram"，则生成直方图并保存为 output/histogram.png。
        """
        # 创建新的 Matplotlib 图形对象，设置图像大小
        # self.config.get("figsize", (10, 6)) 从配置获取图像尺寸，默认值 (10, 6)
        plt.figure(figsize=self.config.get("figsize", (10, 6)))
        
        # data=self.data: 输入的 Pandas DataFrame
        # x="Frequency": 使用 'Frequency' 列作为 X 轴数据（编辑频率，0 到 1）
        # bins=self.config.get("bins", 50): 分箱数量，默认 50，可通过 config 自定义
        sns.histplot(data=self.data, x="Frequency", bins=self.config.get("bins", 50))
        
        # 设置图表标题
        plt.title(self.config.get("title", "RNA editing frequency distribution"))
        plt.xlabel("Frequency")
        plt.ylabel("Count")
        
        self.save(self.config.get("format", "png"))