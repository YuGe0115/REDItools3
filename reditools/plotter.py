# 导入必要的 Python 模块
from abc import ABC, abstractmethod
import pandas as pd
import matplotlib.pyplot as plt

class Plotter(ABC):    
    def __init__(self, data: pd.DataFrame, output_path: str, config: dict):
        """
        初始化 Plotter 类的实例。

        参数:
            data (pd.DataFrame): 包含数据的 Pandas DataFrame，通常是从 REDItools3 的 TSV 文件加载的表格数据。
                                例如，包含 'Region', 'Position', 'Frequency' 等字段。
            output_path (str): 图表输出文件的路径（不含文件扩展名），例如 'output/plot'。
                              最终文件会以指定的格式（如 .png）保存。
            config (dict): 配置字典，包含绘图参数（如图像大小、分辨率、标题等）。
                          示例: {"figsize": [12, 6], "dpi": 300, "title": "RNA 编辑图"}。
                          如果未提供某些配置项，会使用默认值。

        功能:
            - 将输入的 data、output_path 和 config 保存为实例变量，以便子类使用。
            - 确保所有子类都有统一的初始化逻辑。
        """
        self.data = data  # 保存输入的 DataFrame 数据
        self.output_path = output_path  # 保存输出文件路径
        self.config = config  # 保存配置字典

    @abstractmethod
    def plot(self):
        """
        抽象方法，要求子类实现具体的绘图逻辑。

        功能:
            - 由子类（如 HistogramPlotter、ManhattanPlotter）实现，用于生成特定类型的图表。
            - 例如，直方图子类会使用 self.data['Frequency'] 绘制频率分布；
              曼哈顿图子类会使用 self.data['Region'] 和 self.data['Position'] 绘制基因组位置分布。
            - 不需要返回值，图表通过 Matplotlib 的 plt 对象生成。

        注意:
            - 子类必须实现此方法，否则实例化时会抛出 TypeError。
            - 可以使用 self.config 获取用户自定义的绘图参数（如标题、颜色）。
        """
        pass  # 抽象方法不实现具体逻辑，仅定义接口

    def save(self, format: str = "png"):
        """
        保存生成的图表到文件。

        参数:
            format (str): 输出文件的格式，默认为 "png"。支持的格式包括 "png"、"pdf"、"jpg" 等。
                         由 Matplotlib 支持的格式决定。

        功能:
            - 使用 Matplotlib 的 plt.savefig 方法将当前图表保存到文件。
            - 文件名由 self.output_path 和 format 组成，例如 'output/plot.png'。
            - 从 self.config 中获取分辨率（dpi），默认值为 300。
            - 保存后关闭 Matplotlib 的图形对象（plt.close()），避免内存泄漏。

        示例:
            - 如果 self.output_path = 'output/plot'，format = 'png'，config = {'dpi': 300}，
              则保存文件为 'output/plot.png'，分辨率为 300 DPI。
        """
        plt.savefig(f"{self.output_path}.{format}", dpi=self.config.get("dpi", 300))
        # plt.savefig: 保存图表到指定路径
        # f"{self.output_path}.{format}": 使用 f-string 构造文件名，例如 'output/plot.png'
        # dpi=self.config.get("dpi", 300): 从 config 字典获取分辨率，若无则使用默认值 300
        plt.close()  # 关闭当前图形，释放内存，避免多个图表叠加或内存占用过高