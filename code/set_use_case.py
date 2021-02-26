
# import built-in modules
import os

# import third-party modules
import pandas as pd
import wget

# import project modules
from utils import get_logger, load_config, write_config

# Path to config yaml
CONFIG = os.path.join(os.path.dirname(__file__), 'setup_config.yaml')
# Define logger
logger = get_logger(name=__name__)


def main():
    """Main, coordinates step up."""
    # Load config
    config = load_config(CONFIG)

    # Check if analysis directories are present
    check_structure(config["dirs"], config["pipeline_config"])

    # Download fastq files
    ftp_info = config["ftp_info"]
    dirs_dict = config["dirs"]
    raw_data_dir = os.path.join(os.path.join(os.path.dirname(__file__), dirs_dict["root_dir"]),
                                dirs_dict["raw_dir"])
    download_fastq_files(ftp_info,raw_data_dir)


def download_fastq_files(ftp_info_file, raw_data_dir):
    """Method downloads fastq files into raw data directory.

    :param ftp_info_file: Path to file with ftp info to guide download of fastq files. Expects tab
    separated file with the following columns: FASTQ | FASTQ_MD5 | PAIRED_FASTQ | PAIRED_FASTQ_MD5
    | NIST_SAMPLE_NAME.
    :type ftp_info_file: str

    :param raw_data_dir: Path to raw data directory where fastq files should be stored.
    :type raw_data_dir: str
    """
    # Load ftp info
    ftp_path = os.path.join(os.path.dirname(__file__), ftp_info_file)
    ftp_df = pd.read_csv(ftp_path, sep="\t")
    tmp_path = os.path.join(raw_data_dir, '0_blablabla_R1.fastq.gz')
    wget.download(ftp_df.iloc[0,0], out=tmp_path)
    tmp_path = os.path.join(raw_data_dir, '0_blablabla_R2.fastq.gz')
    wget.download(ftp_df.iloc[0,2], out=tmp_path)




def check_structure(dirs, pipe_config):
    """Method makes sure that all required directories are present, namely:
    ngs_mapping, variant_calling, variant_export.

    :param dirs: Dictionary of expected directories associated with the analysis.
    :type dirs: dict

    :param pipe_config: Dictionary with pipeline config yaml.
    :type pipe_config: dict
    """
    root = dirs["root_dir"]
    root_path = os.path.join(os.path.dirname(__file__), root)

    # Raw data directory
    raw_dir = dirs["raw_dir"]
    full_dir_path = os.path.join(root_path, raw_dir)
    if not os.path.exists(full_dir_path):
        logger.warning(f"Directory '{raw_dir}' was not found and was created.")
        # Create directory
        os.makedirs(full_dir_path)

    # Iterate over directories required for analysis
    for directory in dirs["expected_dirs"]:
        full_dir_path = os.path.join(root_path, directory)
        if not os.path.exists(full_dir_path):
            logger.warning(f"Directory '{directory}' was not found and will be created.")
            # Create directory
            os.makedirs(full_dir_path)
            # Update and write config dict
            path = os.path.join(full_dir_path, "config.yaml")
            pipe_config['pipeline_step'].update({'name': directory})
            write_config(path, pipe_config)
    # Log
    logger.info("Required directory and configuration file structured checked.")




if __name__ == "__main__":
    main()
