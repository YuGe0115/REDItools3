import json

def load_config(config_path: str) -> dict:
    """
    加载 JSON 配置文件。

    参数:
        config_path (str): JSON 文件路径。

    返回:
        dict: 配置字典。

    异常:
        ValueError: 如果无法读取配置文件。
    """
    try:
        with open(config_path, "r") as f:
            return json.load(f)
    except Exception as e:
        raise ValueError(f"Unable to read config file: {e}")