
# import built-in modules
import hashlib
import logging

# import third-party modules
import yaml


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


def write_config(path, pipe_config):
    """
    :param path: Path where new config file should be created.
    :type path: str

    :param pipe_config: Dictionary with pipeline configuration.
    :type pipe_config: dict
    """
    with open(path, 'w') as yaml_file:
        yaml.dump(pipe_config, yaml_file, default_flow_style=False)


def load_config(config_path):
    """
    :param config_path: Path to yaml configuration file.
    :type config_path: str

    :return: Returns settings as found in configuration file.
    """
    with open(config_path, "r") as yaml_file:
        cfg = yaml.safe_load(yaml_file)
    return cfg


def file_md5(input_file):
    """
    :param input_file: Path to input file.
    :type input_file: str

    :return: Returns the encoded data in the inputted file in hexadecimal format.
    """
    with open(input_file, 'rb') as f:
        data = f.read()
        return hashlib.md5(data).hexdigest()
