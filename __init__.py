import warnings
warnings.filterwarnings("ignore")

import logging
"""
Configure the logger level for this module.
"""
logging.basicConfig(
    level=logging.ERROR,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

from .utils.temp_manager import TempMangager


temp_manager = TempMangager("./temp")
temp_manager.prepare()  # 创建并清空临时文件夹