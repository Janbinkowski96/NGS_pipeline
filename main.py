# -*- coding: utf-8 -*-S

import click
from source.shell import Shell
from source.file_menager import FileMenager
from utils.utils import save_raw_file
from utils.utils import break_, clear_, print_, help_, check_flags_

@click.command()
@click.option('--help', '-h',
              is_flag=True,
              expose_value=False,
              is_eager=False,
              callback=help_,
              help="Print help message")

@click.option('--read-1', '-1', type=click.Path(exists=True), help='Read 1 path.')
@click.option('--read-2', '-2', type=click.Path(exists=True), help='Read 2 path.')
@click.option('--reference', '-r', type=click.Path(), help='Reference genome index.')
@click.option('--reference-genome', '-g', type=click.Path(), help='Reference genome for VC.')
@click.option('--regions', '-s', type=click.Path(), help='Path to BED file with specyfic regions.')
@click.option('--output', '-o', type=click.Path(), help='Output path')
@click.option('--project-name', '-p', type=str, help='Project name.')
@click.pass_context
def main(ctx, read_1, read_2, reference, reference_genome, regions, output, project_name):
    check_flags_(ctx, read_1, read_2, reference, reference_genome, regions, output, project_name)
    clear_()
    
    main_directory = FileMenager.join_paths(output, project_name)
    status = FileMenager.create_directory(path=main_directory)
    
    if status != 0: # TODO None < 0
        break_("Project directory already exists, to prevent overwrite existing results change project_name.")
    else:
        print_(f"Project directory '{main_directory}' created.", bold=True)
    
    file_menager = FileMenager(main_directory=main_directory)
    shell = Shell(main_directory=main_directory)
    
    #QC dir
    hisat_files = file_menager.join_paths(main_directory, "QC")
    file_menager.create_directory(hisat_files)
    shell.set_path("QC")
    
    # QC
    command = shell.prepare_command_(f"fastp --in1 {read_1} --in2 {read_2}")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    shell.reset_cwd()
    
    # Create aligment dir
    hisat_files = file_menager.join_paths(main_directory, "alignment")
    file_menager.create_directory(hisat_files)
    shell.set_path("alignment")
    
    # Run Hisat2
    command = shell.prepare_command_(f"""hisat2 -x {reference}/ -q 
                                     -1 {read_1} -2 {read_2} 
                                     --summary-file alignment_raport.txt -p 4
                                     --no-spliced-alignment 
                                     -S raw_alignment.sam""")
    print_(command, bold=True, is_command=True)
    print_(str(command).join(""))
    shell.run(command)

    # Convert SAM to BAM
    command = shell.prepare_command_(f"""samtools view raw_alignment.sam -L {regions} -b -S -o raw_alignment.bam""")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    # Delete SAM file
    command = shell.prepare_command_("rm raw_alignment.sam")
    print_(command, bold=True, is_command=True)
    shell.run(command)

    # Samtools fixmate
    command = shell.prepare_command_("samtools fixmate -m raw_alignment.bam fixmate_alignment.bam")
    print_(command, bold=True, is_command=True)
    shell.run(command)

    # Samtools sort
    command = shell.prepare_command_("samtools sort -o sorted_alignment.bam -O BAM fixmate_alignment.bam")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    # Remove duplications
    command = shell.prepare_command_(f"samtools markdup -r -f markduplications.txt sorted_alignment.bam {project_name}_final.bam -s markdup_report.txt")
    print_(command, bold=True, is_command=True)
    shell.run(command)

    # BAM file stats
    command = shell.prepare_command_(f"samtools flagstats -O tsv {project_name}_final.bam")
    print_(command, bold=True, is_command=True)
    stats = shell.run(command)
    save_raw_file(stats, shell.cwd, "final_BAM_stats")
    
    # BAM file indexing
    command = shell.prepare_command_(f"samtools index -b  {project_name}_final.bam")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    # Coverage & depth
    # Create new directory for important BAM stats
    bam_stats_directory = FileMenager.join_paths(main_directory, "depth_coverage_stats")
    FileMenager.create_directory(path=bam_stats_directory)
    shell.reset_cwd()
    shell.set_path(bam_stats_directory)

    # Coverage
    command = shell.prepare_command_(f"bedtools coverage -b ../alignment/{project_name}_final.bam -a {regions}")
    print_(command, bold=True, is_command=True)
    depth_stats = shell.run(command)
    save_raw_file(depth_stats, shell.cwd, "coverage")
    
    # Variant Calling
    # Create new directory for Variant Calling data
    vc_directory = FileMenager.join_paths(main_directory, "VC")
    FileMenager.create_directory(path=vc_directory)
    shell.reset_cwd()
    shell.set_path(vc_directory)
    
    # Mpileup
    command = shell.prepare_command_(f"""bcftools mpileup -E -a DP -a SP -a AD -P ILLUMINA -pm3 -F0.2 -C50 -d 1400000
                                     -f {reference_genome}
                                     -o raw.vcf -O v
                                     ../alignment/{project_name}_final.bam""")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    # Calling
    command = shell.prepare_command_("bcftools call -O v -m -v --ploidy GRCh37 raw.vcf -o variants.vcf")
    print_(command, bold=True, is_command=True)
    shell.run(command)

if __name__ == "__main__":
    main()
