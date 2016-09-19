import logging
from MEOutputFormatter import MEOutputFormatter as mof

class MELoggingFormatter(logging.Formatter):
    """
    This module has been done to modify the way in which the log is written to
    the standard output and the log file.
    """
    FORMATS={\
        "CONFIG":"%(asctime)s - {0}{1}%(levelname)s{2}:\t%(message)s".format(mof.BOLD, mof.BLUE,mof.END),\
        "ERROR": "%(asctime)s - {0}{1}%(levelname)s (%(module)s){2}{0}:\t%(message)s{2}".format(mof.BOLD, mof.RED,mof.END),\
        "WARNING":"%(asctime)s - {0}{1}%(levelname)s{2}:\t%(message)s".format(mof.BOLD, mof.YELLOW,mof.END),\
        "INFO": "%(asctime)s - {0}{1}%(levelname)s{2}:\t%(message)s".format(mof.BOLD, mof.GREEN,mof.END),\
        "DEBUG":"%(asctime)s - {0}{1}%(levelname)s{2} (%(module)s|%(funcName)s:%(lineno)d):\t%(message)s".format(mof.BOLD, mof.PURPLE,mof.END),\
        "DEFAULT":"%(asctime)s - {0}%(levelname)s{1}:\t%(message)s".format(mof.BOLD,mof.END)\
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
