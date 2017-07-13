import logging

LOG_LEVEL_CHOICES=["DEBUG","INFO","WARNING","ERROR"]
CONFIG_LEVEL=45

class MEOutputFormatter:
    """
    This module incorporates different constant color values from bash shell
    for a "pretty" visualization of the log.
    """
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

class MELoggingFormatter(logging.Formatter):
    """
    This module has been done to modify the way in which the log is written to
    the standard output and the log file.
    """
    FORMATS={\
        "CONFIG":"%(asctime)s - {0}{1}%(levelname)s{2}:\t%(message)s".format(\
            MEOutputFormatter.BOLD, MEOutputFormatter.BLUE,MEOutputFormatter.END),\
        "ERROR": "%(asctime)s - {0}{1}%(levelname)s (%(module)s){2}{0}:\t%(message)s{2}".format(\
            MEOutputFormatter.BOLD, MEOutputFormatter.RED,MEOutputFormatter.END),\
        "WARNING":"%(asctime)s - {0}{1}%(levelname)s{2}:\t%(message)s".format(\
            MEOutputFormatter.BOLD, MEOutputFormatter.YELLOW,MEOutputFormatter.END),\
        "INFO": "%(asctime)s - {0}{1}%(levelname)s{2}:\t%(message)s".format(\
            MEOutputFormatter.BOLD, MEOutputFormatter.GREEN,MEOutputFormatter.END),\
        "DEBUG":"%(asctime)s - {0}{1}%(levelname)s{2} (%(module)s|%(funcName)s:%(lineno)d):\t%(message)s".format(\
            MEOutputFormatter.BOLD, MEOutputFormatter.PURPLE,MEOutputFormatter.END),\
        "DEFAULT":"%(asctime)s - {0}%(levelname)s{1}:\t%(message)s".format(\
            MEOutputFormatter.BOLD,MEOutputFormatter.END)\
       }

    def __init__(self,fmt,datefmt):
        logging.Formatter.__init__(self, fmt, datefmt)

    def format(self, record):
        original_fmt=self._fmt
        try:
            self._fmt = MELoggingFormatter.FORMATS[record.levelname]
        except:
            self._fmt = MELoggingFormatter.FORMATS["DEFAULT"]

        return logging.Formatter.format(self, record)
