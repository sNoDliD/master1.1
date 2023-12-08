import logging
import pathlib
import sys
from logging.handlers import RotatingFileHandler

logger = logging.getLogger("bot.logger")
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter(
    '%(levelname)-6s [%(asctime)s] %(message)-100s      <%(filename)s:%(lineno)d ~%(threadName)s~>')

path = pathlib.Path(__file__).parent.resolve() / "bot.log"
fh = RotatingFileHandler(path, maxBytes=30 * 1024 * 1024, backupCount=1, encoding="utf-8")
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler(stream=sys.stdout)
ch.setFormatter(formatter)
logger.addHandler(ch)
