#!/usr/bin/env python3
"""
Direct bam2seqz wrapper — calls sequenza Python internals to bypass
the Python 3.13 argparse CLI conflict in sequenza-utils.
"""
import sys, os, gzip, argparse
sys.path.insert(0, '/home/ec2-user/miniforge3/lib/python3.13/site-packages')

from sequenza.programs.bam2seqz import bam2seqz
from sequenza.misc import SeqzLogger

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--normal',  required=True)
    p.add_argument('--tumor',   required=True)
    p.add_argument('--fasta',   required=True)
    p.add_argument('--gc_file', required=True)
    p.add_argument('--output',  required=True)
    p.add_argument('--parallel', type=int, default=8)
    p.add_argument('--chromosome', default=None)
    args = p.parse_args()

    extra = [
        '-n', args.normal,
        '-t', args.tumor,
        '-F', args.fasta,
        '--gc_file', args.gc_file,
        '-o', args.output,
        '--parallel', str(args.parallel),
    ]
    if args.chromosome:
        extra += ['--chromosome', args.chromosome]

    import argparse as _ap
    subparsers = _ap.ArgumentParser().add_subparsers()
    log = SeqzLogger(level=30)
    bam2seqz(subparsers, 'bam2seqz', extra, log)

if __name__ == '__main__':
    main()
