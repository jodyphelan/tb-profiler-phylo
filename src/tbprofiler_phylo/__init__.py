"""
Module for tb-profiler plugin that
provides phylogenetic functionality.
"""

__version__ = "0.0.1"

from tbprofiler import ProfilePlugin
from tbprofiler.consensus import get_consensus_vcf, prepare_sample_consensus
from pathogenprofiler.utils import run_cmd
import argparse
import logging
import os
import filelock
from joblib import Parallel,delayed
from tqdm import tqdm





def usher_add_sample(args: argparse.Namespace) -> None:
    logging.info("Adding sample to phylogeny")


    if args.vcf:
        args.wg_vcf = args.vcf
    else:
        args.wg_vcf = args.files_prefix + ".vcf.gz"

    args.tmp_masked_vcf = f"{args.files_prefix}.masked.vcf.gz"
    args.input_phylo = f"{args.dir}/results/phylo.pb"
    args.tmp_output_phylo = f"{args.files_prefix}.pb"
    args.output_nwk = f"{args.files_prefix}.nwk"

    if not os.path.isfile(args.input_phylo):
        logging.error("Phylogeny doesn't exist. Please create one first with `tb-profiler-tools`")
        quit("Exiting!")


    lock = filelock.SoftFileLock(args.input_phylo + ".lock")

    cwd = os.getcwd()
    args.tmp_masked_vcf = get_consensus_vcf(args.prefix, args.wg_vcf,args)
    with lock:
        os.chdir(args.temp)

        run_cmd("usher --vcf %(tmp_masked_vcf)s --load-mutation-annotated-tree %(input_phylo)s --save-mutation-annotated-tree %(tmp_output_phylo)s --write-uncondensed-final-tree" % vars(args))
        run_cmd("mv uncondensed-final-tree.nh %(output_nwk)s" % vars(args))
        for f in ["mutation-paths.txt","placement_stats.tsv"]:
            if os.path.exists(f):
                os.remove(f)
        run_cmd("mv %(tmp_output_phylo)s %(input_phylo)s " % vars(args))
        os.chdir(cwd)



def prepare_usher(treefile: str,vcf_file: str) -> None:
    run_cmd(f"usher --tree {treefile} --vcf {vcf_file} --collapse-tree --save-mutation-annotated-tree phylo.pb")



def wrapper_function(s: str,args: argparse.Namespace) -> str:
    args.bam = f"{args.dir}/bam/{s}.bam"
    return prepare_sample_consensus(
        sample_name=s,
        ref=args.conf['ref'],
        input_vcf=f"{args.dir}/vcf/{s}.vcf.gz",
        output_file=f"{args.temp}/{s}.consensus.fasta",
        excluded_regions=args.conf['bedmask']
    )

def calculate_phylogeny(args: argparse.Namespace) -> None:
    samples = [l.strip() for l in open(args.samples)]
    args.tmp_masked_vcf = f"{args.files_prefix}.masked.vcf.gz"
    
    alignment_file = f"{args.files_prefix}.aln"
    consensus_files = [r for r in tqdm(Parallel(n_jobs=args.threads,return_as='generator')(delayed(wrapper_function)(s,args) for s in samples),desc="Generating consensus sequences",total=len(samples))]

    run_cmd(f"cat {' '.join(consensus_files)} > {alignment_file}")
    alignment_file_plus_ref = f"{args.files_prefix}.aln.plus_ref"
    run_cmd(f"cat {args.conf['ref']} > {alignment_file_plus_ref}")
    run_cmd(f"cat {alignment_file} >> {alignment_file_plus_ref}")
    tmp_vcf = f"{args.files_prefix}.vcf"
    run_cmd(f"faToVcf {alignment_file_plus_ref} {tmp_vcf}")
    run_cmd(f"iqtree -s {alignment_file} -m GTR+G -nt 2",desc="Running IQTree")
    prepare_usher(f"{alignment_file}.treefile",tmp_vcf)
    run_cmd(f"mv phylo.pb {args.dir}/results/")
    os.remove("condensed-tree.nh")

class PhyloPlugin(ProfilePlugin):
    __cli_params__ = [
        {
            "args": ["--update-phylo"],
            'kwargs': {
                "action": "store_true",
                "help": "Add sample to phylogeny"
            }
        }
    ]
    def run(self,args, _ ) -> None:
        if args.update_phylo:
            usher_add_sample(args)
        else:
            pass
            
    def pre_process(self,args) -> None:
        if args.update_phylo:
            logging.debug("Setting `--call-whole-genome` to True")
            args.call_whole_genome = True

    def post_process(self,args) -> None:
        pass