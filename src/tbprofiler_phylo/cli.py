import tbprofiler as tbp
import argparse
import sys
import os
from uuid import uuid4
import pathogenprofiler as pp
from glob import glob
import logging
from rich.logging import RichHandler
import tbprofiler_phylo

__softwarename__ = 'tbprofiler'

def cleanup(files_prefix):
    for f in glob(f"{files_prefix}*"):
        os.remove(f)

def main_phylogeny(args):
    args.conf['variant_filters'] = {
        'depth_soft': args.min_dp
    }
    tbprofiler_phylo.calculate_phylogeny(args)
    cleanup(args.files_prefix)

def prepare_usher(args):
    tbprofiler_phylo.prepare_usher(args.tree,args.vcf)

def cli():
    argparser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = argparser.add_subparsers(help='commands')
    parser_sub = subparsers.add_parser('create_phylogeny', help='Calculate phylogeny')
    parser_sub.add_argument('--samples',required=True,help="Samples files")
    parser_sub.add_argument('--temp',help="Temp firectory to process all files",type=str,default=".")
    parser_sub.add_argument('--db',default='tbdb',help='Mutation panel name')
    parser_sub.add_argument('--external-db',type=str,help='Path to db files prefix (overrides "--db" parameter)')
    parser_sub.add_argument('--dir','-d',default=".",help='Storage directory')
    parser_sub.add_argument('--db-dir',default=tbp.utils.get_default_db_dir(),help='DB directory')
    parser_sub.add_argument('--min-dp',default=10,type=int,help='Bases with depth below this cutoff will be marked as missing')
    parser_sub.add_argument('--threads',default=1,type=int,help='Total number of threads')
    parser_sub.add_argument('--logging',default="INFO",choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"],help='Logging level')
    parser_sub.set_defaults(func=main_phylogeny)

    parser_sub = subparsers.add_parser('prepare_usher', help='Update tbdb')
    parser_sub.add_argument('--tree',required=True,help="Samples files")
    parser_sub.add_argument('--vcf',help="Temp firectory to process all files",required=True)
    parser_sub.add_argument('--temp',help="Temp firectory to process all files",type=str,default=".")
    parser_sub.add_argument('--dir','-d',default=".",help='Storage directory')
    parser_sub.add_argument('--db-dir',default=tbp.utils.get_default_db_dir(),help='DB directory')
    parser_sub.add_argument('--logging',default="INFO",choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"],help='Logging level')
    parser_sub.set_defaults(func=prepare_usher)

    args = argparser.parse_args()
    if hasattr(args, 'func'):

        logging.basicConfig(
            level=args.logging, format="%(message)s", datefmt="[%X]", handlers=[RichHandler()]
        )

        args.software_name = __softwarename__
        args.tmp_prefix = str(uuid4())
        args.files_prefix = os.path.abspath(f"{args.temp}/{args.tmp_prefix}")
        if hasattr(args,'dir'):
            args.dir = os.path.abspath(args.dir)
        if hasattr(args, 'db'):
            if args.db=="who-v2" and not args.external_db and pp.nofile(sys.base_prefix+"/share/tbprofiler/who-v2.fasta"):
                logging.error("Can't find the tbdb file at %s. Please run 'tb-profiler update_tbdb' to load the default library or specify another using the '--external_db' flag" % sys.base_prefix)
                quit(1)
            args.conf = pp.get_db(tbp.utils.get_default_data_dir(),args.external_db)
        args.func(args)
    else:
        argparser.print_help(sys.stderr)