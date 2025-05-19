#!/usr/bin/env python3
"""
Filter_TTAA.py
-------------
Filter FASTQ reads to keep those with TTAA motifs, which is 
expected after adapter trimming in PiggyBac transposon experiments.

This enhanced version checks for:
1. TTAA at the beginning of the read (5' end, forward strand)
2. TTAA at the end of the read (3' end, forward strand)
3. AATT at the end of the read (3' end, reverse complement)

This script also filters reads by minimum length, which is useful
after trimming as some reads may become too short.

Usage Example:
  python Filter_TTAA.py \
    --fastq bc07_bc01_all_trimmed.fq \
    --min_length 50 \
    --out bc07_bc01_all_trimmed_TTAA.fq
"""
import sys
import re
import argparse
import gzip
import time
from datetime import datetime

def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTQ reads to keep only those containing TTAA motifs."
    )
    parser.add_argument("--fastq", required=True,
        help="Input FASTQ file (can be .gz or uncompressed).")
    parser.add_argument("--out", required=True,
        help="Output FASTQ (will be uncompressed; user can gzip later).")
    parser.add_argument("--min_length", type=int, default=50,
        help="Minimum read length to keep (default: 50).")
    parser.add_argument("--report", default=None,
        help="Optional path to write detailed filtering report.")
    parser.add_argument("--skip_end_checks", action="store_true",
        help="Skip checking for TTAA/AATT at read ends (only check 5' start)")

    args = parser.parse_args()
    filter_start_time = time.time()
    
    # Open input file
    try:
        if args.fastq.endswith(".gz"):
            fastq_in = gzip.open(args.fastq, "rt")
        else:
            fastq_in = open(args.fastq, "r")
    except IOError as e:
        log(f"ERROR opening input FASTQ: {str(e)}")
        sys.exit(1)

    # Open output file
    try:
        fastq_out = open(args.out, "w")
    except IOError as e:
        log(f"ERROR opening output file: {str(e)}")
        fastq_in.close()
        sys.exit(1)

    log(f"Filtering reads from {args.fastq} -> {args.out}")
    log(f"Keeping reads with TTAA motifs and length >= {args.min_length}")
    
    check_both_ends = not args.skip_end_checks
    if check_both_ends:
        log("Checking for TTAA at both 5' and 3' ends of reads")
    else:
        log("Only checking for TTAA at 5' end of reads")

    # Statistics tracking
    stats = {
        "total_reads": 0,
        "passed_ttaa_start": 0,  # TTAA at start (5' end)
        "passed_ttaa_end": 0,    # TTAA at end (3' end)
        "passed_aatt_end": 0,    # AATT at end (3' end, reverse complement)
        "failed_ttaa": 0,        # No TTAA found at either end
        "too_short": 0           # Read length too short
    }

    # Process FASTQ in chunks of 4 lines for better performance
    try:
        chunk = []
        for line in fastq_in:
            chunk.append(line)
            
            if len(chunk) == 4:
                stats["total_reads"] += 1
                
                # Progress indicator for every 100k reads
                if stats["total_reads"] % 100000 == 0:
                    elapsed = time.time() - filter_start_time
                    log(f"Processed {stats['total_reads']:,} reads ({elapsed:.1f} seconds)")
                
                # Process the chunk (4 lines = 1 read)
                header_line = chunk[0]
                seq_line = chunk[1]
                plus_line = chunk[2]
                qual_line = chunk[3]
                
                # Get the sequence (removing newline)
                seq = seq_line.strip()
                
                # Check sequence length
                if len(seq) < args.min_length:
                    stats["too_short"] += 1
                    chunk = []
                    continue
                
                # Check for TTAA/AATT motifs
                has_ttaa = False
                
                # Check start of sequence for TTAA
                if seq.startswith("TTAA"):
                    stats["passed_ttaa_start"] += 1
                    has_ttaa = True
                
                # Check end of sequence if enabled
                if check_both_ends:
                    # Check for TTAA at end
                    if seq.endswith("TTAA"):
                        stats["passed_ttaa_end"] += 1
                        has_ttaa = True
                    
                    # Check for AATT at end (reverse complement of TTAA)
                    elif seq.endswith("AATT"):
                        stats["passed_aatt_end"] += 1
                        has_ttaa = True
                
                # Write read to output if it has TTAA/AATT
                if has_ttaa:
                    fastq_out.write("".join(chunk))
                else:
                    stats["failed_ttaa"] += 1
                
                # Reset for next chunk
                chunk = []
    except Exception as e:
        log(f"ERROR during FASTQ processing: {str(e)}")
    finally:
        fastq_in.close()
        fastq_out.close()

    # Report statistics
    elapsed = time.time() - filter_start_time
    log(f"Filtering complete in {elapsed:.1f} seconds")
    log(f"Total reads processed: {stats['total_reads']:,}")
    
    total_examined = stats["total_reads"] - stats["too_short"]
    if total_examined > 0:
        # Simpler calculation of total passing reads
        total_passed = total_examined - stats["failed_ttaa"]
        total_pass_percent = (total_passed / total_examined) * 100
        
        log(f"Reads passing any TTAA filter: {total_passed:,} ({total_pass_percent:.1f}%)")
        log(f"  - TTAA at 5' end: {stats['passed_ttaa_start']:,} ({stats['passed_ttaa_start']/total_examined*100:.1f}%)")
        
        if check_both_ends:
            log(f"  - TTAA at 3' end: {stats['passed_ttaa_end']:,} ({stats['passed_ttaa_end']/total_examined*100:.1f}%)")
            log(f"  - AATT at 3' end: {stats['passed_aatt_end']:,} ({stats['passed_aatt_end']/total_examined*100:.1f}%)")
        
        log(f"Reads failing TTAA filter: {stats['failed_ttaa']:,} ({stats['failed_ttaa']/total_examined*100:.1f}%)")
    
    log(f"Reads too short (<{args.min_length}bp): {stats['too_short']:,} ({stats['too_short']/stats['total_reads']*100:.1f}%)")
    
    # Write detailed report if requested
    if args.report:
        try:
            with open(args.report, "w") as report_file:
                report_file.write("# TTAA Filtering Report\n")
                report_file.write(f"Input file: {args.fastq}\n")
                report_file.write(f"Output file: {args.out}\n")
                report_file.write(f"Minimum length: {args.min_length}\n")
                report_file.write(f"Check both ends: {check_both_ends}\n")
                report_file.write(f"Date: {datetime.now()}\n\n")
                
                report_file.write("## Statistics\n")
                report_file.write(f"Total reads processed: {stats['total_reads']:,}\n")
                
                total_examined = stats["total_reads"] - stats["too_short"]
                if total_examined > 0:
                    # Use the same simplified calculation in the report
                    total_passed = total_examined - stats["failed_ttaa"] 
                    total_pass_percent = (total_passed / total_examined) * 100
                    
                    report_file.write(f"Reads passing any TTAA filter: {total_passed:,} ({total_pass_percent:.1f}%)\n")
                    report_file.write(f"  - TTAA at 5' end: {stats['passed_ttaa_start']:,} ({stats['passed_ttaa_start']/total_examined*100:.1f}%)\n")
                    
                    if check_both_ends:
                        report_file.write(f"  - TTAA at 3' end: {stats['passed_ttaa_end']:,} ({stats['passed_ttaa_end']/total_examined*100:.1f}%)\n")
                        report_file.write(f"  - AATT at 3' end: {stats['passed_aatt_end']:,} ({stats['passed_aatt_end']/total_examined*100:.1f}%)\n")
                    
                    report_file.write(f"Reads failing TTAA filter: {stats['failed_ttaa']:,} ({stats['failed_ttaa']/total_examined*100:.1f}%)\n")
                
                report_file.write(f"Reads too short (<{args.min_length}bp): {stats['too_short']:,} ({stats['too_short']/stats['total_reads']*100:.1f}%)\n")
            
            log(f"Detailed report written to {args.report}")
        except IOError as e:
            log(f"ERROR writing report: {str(e)}")

if __name__ == "__main__":
    main()
