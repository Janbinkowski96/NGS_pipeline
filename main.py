# -*- coding: utf-8 -*-

import click
from utils.utils import check_directory_
from utils.utils import update_read_path_
from utils.utils import create_results_dir_
from utils.utils import break_, clear_, print_
from utils.source import alignment, sam_to_bam

@click.command()
@click.option('--read-1', type=click.Path(exists=True), help='Read 1 path.', default="resources/519860_S6_L001_R1_001.fastq")
@click.option('--read-2', type=click.Path(exists=True), help='Read 2 path.', default="resources/519860_S6_L001_R2_001.fastq")
@click.option('--reference', type=click.Path(), help='Reference.', default="resources/hisat-index-builded/")
@click.option('--project-name', type=str, help='Project name.', default="TEST")

@click.option('--gunzip', '-g', is_flag=True, help='Use if reads are gzpi-ed.', default=False)

def main(read_1, read_2, reference, gunzip, project_name):
    raw_elements_dir, final_project_dir = create_results_dir_(project_name)

    if gunzip:
        gunzip_(read_1)
        gunzip_(read_2)
        read_1 = update_read_path_(read_1)
        read_2 = update_read_path_(read_2)
    
    sam_file_path = alignment(reference, read_1, read_2, raw_elements_dir)
    sam_to_bam(raw_elements_dir, sam_file_path)
    
if __name__ == "__main__":
    main()