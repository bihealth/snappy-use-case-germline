
# import built-in modules
import logging


def get_logger(name):
    """
    :return: Returns logger.
    """
    logger = logging.getLogger(name=name)
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    format_str = '%(asctime)s - %(filename)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s'
    formatter = logging.Formatter(format_str)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Return
    return logger
