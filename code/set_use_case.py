
# import built-in modules
import csv
import gzip
import os
import re

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
    dirs_dict = config["dirs"]
    root_dir_path = os.path.join(os.path.dirname(__file__), dirs_dict["root_dir"])

    # Check if analysis directories are present
    check_structure(config["dirs"], config["pipeline_config"])
    # Path to bed files
    bed_file = os.path.join(root_dir_path, dirs_dict["bed_file"])
    # Path fastq files
    raw_data_dir = os.path.abspath(os.path.join(root_dir_path, dirs_dict["raw_dir"]))

    # Update analysis config
    analysis_config = os.path.join(root_dir_path, config["analysis_config"])
    update_analysis_config(analysis_config, raw_data_dir, bed_file)

def update_analysis_config(analysis_config, raw_data_dir, bed_file):
    """Method updates analysis config with missing paths.

    :param analysis_config: Path configuration file used in the analysis.
    :type analysis_config: str

    :param raw_data_dir: Full path to raw data directory.
    :type raw_data_dir: str

    :param bed_file: Full path to kit bed file.
    :type bed_file: str
    """
    logger.info(f"Replacing entries in analysis config file:\n"
                f"- bed file:{bed_file}\n"
                f"- raw data dir: {raw_data_dir}")
    # Load config
    config = load_config(analysis_config)

    # Modify bed file path
    target_coverage_report_dict = config['step_config']['ngs_mapping']['target_coverage_report']
    target_coverage_report_dict['path_target_interval_list_mapping'][0]['path'] = bed_file
    config['step_config']['ngs_mapping']['target_coverage_report'] = target_coverage_report_dict

    # Modify raw data dir path
    config['data_sets']['mundlos_limb']['search_paths'] = [raw_data_dir]

    # Replace
    write_config(analysis_config, config)


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

    # Raw and static data directory
    data_dir = [os.path.join(root_path, dirs["raw_dir"]),
                os.path.join(root_path, dirs["static_dir"])]
    # Iterate
    for full_dir_path in data_dir:
        if not os.path.exists(full_dir_path):
            dir_name = os.path.basename(full_dir_path)
            logger.warning(f"Directory '{dir_name}' was not found and was created.")
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
