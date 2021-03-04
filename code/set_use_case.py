
# import built-in modules
import csv
import os

# import third-party modules
import wget

# import project modules
from utils import get_logger, load_config, file_md5, write_config

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
    coordinate_fastq_download(ftp_info, raw_data_dir)


def coordinate_fastq_download(ftp_info_file, raw_data_dir):
    """Method coordinates the fastq files download into raw data directory.

    :param ftp_info_file: Path to file with ftp info to guide download of fastq files. Expects tab
    separated file with the following columns: FASTQ | FASTQ_MD5 | PAIRED_FASTQ | PAIRED_FASTQ_MD5
    | NIST_SAMPLE_NAME.
    :type ftp_info_file: str

    :param raw_data_dir: Path to raw data directory where fastq files should be stored.
    :type raw_data_dir: str
    """
    # Initialise variables
    r1_ftp_url = None
    r1_ftp_md5 = None
    r1_name = None
    r2_ftp_url = None
    r2_ftp_md5 = None
    r2_name = None

    # Load ftp info
    ftp_path = os.path.join(os.path.dirname(__file__), ftp_info_file)
    with open(ftp_path) as tsv:
        for line in csv.reader(tsv, delimiter="\t"):
            # Ignores header
            if 'FASTQ' in line:
                continue
            # Only processes the first fastq file; break after
            else:
                # R1 info
                r1_ftp_url = line[0]
                r1_ftp_md5 = line[1]
                r1_name = r1_ftp_url.split('/')[-1]
                # R2 info
                r2_ftp_url = line[2]
                r2_ftp_md5 = line[3]
                r2_name = r2_ftp_url.split('/')[-1]
                break

    # Download R1
    download_and_validate(ftp_url=r1_ftp_url, file_name=r1_name,
                          file_expected_md5=r1_ftp_md5, out_dir=raw_data_dir)
    # Download R2
    download_and_validate(ftp_url=r2_ftp_url, file_name=r2_name,
                          file_expected_md5=r2_ftp_md5, out_dir=raw_data_dir)


def download_and_validate(ftp_url, file_name, file_expected_md5, out_dir):
    """Method downloads file and validates the download using the expected md5.
    If expected and observed md5 values are different, it will raise Exception.

    :param ftp_url: FTP url to file.
    :type ftp_url: str

    :param file_name: File name, will be used to name file locally.
    :type file_name: str

    :param file_expected_md5: File expected md5, will be used to validate download.

    :param out_dir: Path to output directory, where file should be saved.
    :type out_dir: str
    """
    # Log
    logger.info(f"Downloading {file_name}...")
    # Download
    local_path = os.path.join(out_dir, file_name)
    if os.path.exists(local_path):
         os.remove(local_path)  # if exist, remove it directly
    wget.download(ftp_url, out=local_path)
    # MD5 file
    file_observed_md5 = file_md5(local_path)

    if file_expected_md5 != file_observed_md5:
        logger.error(f"File: {file_name}; ftp md5 {file_expected_md5}; "
                     f"local md5 {file_observed_md5}.")
        raise Exception("Download failed: md5 observed different from expected.")

    # Log
    logger.info("done!")


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
