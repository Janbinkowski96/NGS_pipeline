# -*- coding: utf-8 -*-S

import click
from source.shell import Shell
from source.file_menager import FileMenager
from utils.utils import save_raw_file
from utils.utils import break_, clear_, print_, help_, check_flags_

from source.config import *

@click.command()
@click.option('--help', '-h',
              is_flag=True,
              expose_value=False,
              is_eager=False,
              callback=help_,
              help="Print help message")

@click.option('--read-1', '-1', type=click.Path(exists=True), help='Read 1 path.')
@click.option('--read-2', '-2', type=click.Path(exists=True), help='Read 2 path.')
@click.option('--index', '-i', type=click.Path(), help='Reference genome index.')
@click.option('--reference-genome', '-g', type=click.Path(), help='Reference genome for VC.')
@click.option('--regions', '-s', type=click.Path(), help='Path to BED file with specyfic regions.')
@click.option('--algorithm', '-a', type=click.Choice(["hisa2", "bowtie2", "bwa"]), help='Alignment algorithm.')
@click.option('--output', '-o', type=click.Path(), help='Output path')
@click.option('--project-name', '-p', type=str, help='Project name.')
@click.pass_context
def main(ctx, read_1, read_2, index, reference_genome, regions, algorithm, output, project_name):
    check_flags_(ctx, read_1, read_2, index, reference_genome, regions, algorithm, output, project_name)
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
    command = shell.prepare_command_(f"fastp -w {THREADS_NUMBER} --in1 {read_1} --in2 {read_2} -o read_1.fastq -O read_2.fastq")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    shell.reset_cwd()

    # Create aligment dir
    alignment_dir = file_menager.join_paths(main_directory, "alignment")
    file_menager.create_directory(alignment_dir)
    shell.set_path("alignment")
    
    # Run Mapper
    if algorithm == "hisat2":
        command = shell.prepare_command_(f"""hisat2 -x {index}/ -q 
                                         -1 ../QC/read_1.fastq -2  ../QC/read_2.fastq
                                         --summary-file alignment_raport.txt -p {THREADS_NUMBER}
                                         --no-spliced-alignment 
                                         -S raw_alignment.sam""")
    elif algorithm == "bowtie2":
        command = shell.prepare_command_(f"""bowtie2 -x {index}/ -1 ../QC/read_1.fastq
                                         -2 ../QC/read_2.fastq --phred33 --sensitive 
                                         -p {THREADS_NUMBER} -S raw_alignment.sam
                                         """)
    else:
    	command = shell.prepare_command_(f"""bwa mem -M -t {THREADS_NUMBER} -Y -o raw_alignment.sam {index}  
    					  ../QC/read_1.fastq ../QC/read_2.fastq
                                         """)
                                         
    print_(command, bold=True, is_command=True)
    print_(str(command).join(""))
    shell.run(command)

    # Convert SAM to BAM
    command = shell.prepare_command_(f"""samtools view raw_alignment.sam -L {regions} -b -S -@ {THREADS_NUMBER} -o raw_alignment.bam""")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    # Delete SAM file
    command = shell.prepare_command_("rm raw_alignment.sam")
    print_(command, bold=True, is_command=True)
    shell.run(command)

    # Samtools fixmate
    command = shell.prepare_command_("samtools fixmate -m raw_alignment.bam fixmate_alignment.bam -@ {THREADS_NUMBER}")
    print_(command, bold=True, is_command=True)
    shell.run(command)

    # Samtools sort
    command = shell.prepare_command_("samtools sort -@ {THREADS_NUMBER} -o sorted_alignment.bam -O BAM fixmate_alignment.bam")
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
    shell.set_path("depth_coverage_stats")

    # Coverage
    command = shell.prepare_command_(f"bedtools coverage -a {regions} -b ../alignment/{project_name}_final.bam")
    print_(command, bold=True, is_command=True)
    depth_stats = shell.run(command)
    save_raw_file(depth_stats, shell.cwd, "coverage")
    
    # Variant Calling
    # Create new directory for Variant Calling data
    vc_directory = FileMenager.join_paths(main_directory, "VC")
    FileMenager.create_directory(path=vc_directory)
    shell.reset_cwd()
    shell.set_path("VC")
    
    # Mpileup
    command = shell.prepare_command_(f"""bcftools mpileup -a DP,SP,AD,ADF,ADR  -d 100 -p -m 2 -F 0.2 -A -C 50
                                     -f {reference_genome}
                                     -o raw.vcf -O v --threads {THREADS_NUMBER}
                                     ../alignment/{project_name}_final.bam""")
                                     
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    # Calling
    command = shell.prepare_command_("bcftools call -cv --threads {THREADS_NUMBER} --ploidy GRCh37 -p 0.09 -Ov -o variants.vcf raw.vcf")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    
    # GZIP and Index VCF
    command = shell.prepare_command_("bgzip variants.vcf")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    command = shell.prepare_command_("bcftools index -t variants.vcf.gz")
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    # Normalization
    command = shell.prepare_command_(f"""bcftools norm variants.vcf.gz -Ov --check-ref w -f {reference_genome} 
    					  -R {regions} -o final.vcf""")
    
    print_(command, bold=True, is_command=True)
    shell.run(command)
    
    # Annotations
    command = shell.prepare_command_(f"""
                                     vep -i variants.vcf -o variants_annotated.json --cache 
                                     -port 3337 -force_overwrite 
                                     --check_existing --symbol -vcf -biotype 
                                     --merged --sift b --everything --json
                                     -fields Uploaded_variation, Location, Allele, 
                                     VARIANT_CLASS,BIOTYPE, Consequence, 
                                     Existing_variation, CLIN_SIG,SYMBOL, 
                                     SYMBOL_SOURCE, Gene, Feature_type, Feature, 
                                     cDNA_position, CDS_position, Protein_position, Amino_acids, 
                                     Codons, IMPACT, DISTANCE, EXON, INTRON, 
                                     DOMAINS, SIFT, PolyPhen, MAX_AF, MAX_AF_POPS,
                                     AF, AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, 
                                     AA_AF, EA_AF, gnomAD_AF, gnomAD_AFR_AF, gnomAD_AMR_AF, 
                                     gnomAD_ASJ_AF, gnomAD_EAS_AF, gnomAD_FIN_AF, 
                                     gnomAD_NFE_AF, gnomAD_OTH_AF, gnomAD_SAS_AF
                                     """)
    #print_(command, bold=True, is_command=True)
    #shell.run(command)
    
if __name__ == "__main__":
    main()
