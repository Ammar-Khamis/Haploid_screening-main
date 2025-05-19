#!/usr/bin/env python3
"""
Improved_trim_fastq.py
----------------------
An enhanced version of the trimming script for alignment-based trimming of FASTQ reads.

Features added:
- Better error handling and reporting
- Progress indication
- Performance optimization
- Statistics on trimming results

Usage Example:
  python Improved_trim_fastq.py \
    --fastq bc07_bc01_all.fq.gz \
    --trimdata bc07_bc01_for_trimdata.tsv \
    --out bc07_bc01_all_trimmed.fq
"""
import sys
import re
import argparse
import gzip
import time
from datetime import datetime
from collections import Counter

def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}")

def main():
    parser = argparse.ArgumentParser(
        description="Trim reads by a custom offset from a TSV (start+length, strand)."
    )
    parser.add_argument("--fastq", required=True,
        help="Input FASTQ file (can be .gz or uncompressed).")
    parser.add_argument("--trimdata", required=True,
        help="TSV with columns for readID, start, length, and strand.")
    parser.add_argument("--out", required=True,
        help="Output FASTQ (will be uncompressed; user can gzip later).")
    parser.add_argument("--report", default=None,
        help="Optional path to write detailed trimming report.")

    args = parser.parse_args()
    trim_start_time = time.time()
    
    # 1) Build the trim dictionary from the TSV
    log(f"Loading trim data from {args.trimdata}...")
    trimdict = {}
    trim_data_stats = {"total_lines": 0, "skipped_lines": 0}
    
    try:
        with open(args.trimdata, "r") as f:
            for line in f:
                trim_data_stats["total_lines"] += 1
                fields = line.strip().split("\t")
                if len(fields) < 6:
                    trim_data_stats["skipped_lines"] += 1
                    continue

                read_id = fields[5].strip()
                start = int(fields[2])
                length = int(fields[3])
                strand = fields[4]

                if strand == "+":
                    v = start + length
                else:
                    v = -1 * (start + length)

                trimdict[read_id] = v
    except (IOError, ValueError) as e:
        log(f"ERROR loading trim data: {str(e)}")
        sys.exit(1)

    log(f"Loaded trim data for {len(trimdict)} read IDs from {args.trimdata}.")
    log(f"Processed {trim_data_stats['total_lines']} lines, skipped {trim_data_stats['skipped_lines']} incomplete lines.")

    # 2) Decide how to open the input FASTQ
    try:
        if args.fastq.endswith(".gz"):
            fastq_in = gzip.open(args.fastq, "rt")
        else:
            fastq_in = open(args.fastq, "r")
    except IOError as e:
        log(f"ERROR opening input FASTQ: {str(e)}")
        sys.exit(1)

    # 3) Open the output uncompressed
    try:
        fastq_out = open(args.out, "w")
    except IOError as e:
        log(f"ERROR opening output file: {str(e)}")
        fastq_in.close()
        sys.exit(1)

    log(f"Trimming reads from {args.fastq} -> {args.out}")

    # Regex to match read IDs: e.g. @3ef3c0ba-...
    fqIDre = re.compile(r"^@(\w{8}-\w{4}-\w{4}-\w{4}-\w{12})")

    # Statistics tracking
    stats = {
        "total_reads": 0,
        "trimmed_left": 0,
        "trimmed_right": 0,
        "not_trimmed": 0,
        "no_match_id": 0,
        "trim_amounts": Counter()
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
                    elapsed = time.time() - trim_start_time
                    log(f"Processed {stats['total_reads']:,} reads ({elapsed:.1f} seconds)")
                
                # Process the chunk (4 lines = 1 read)
                header_line = chunk[0]
                seq_line = chunk[1]
                plus_line = chunk[2]
                qual_line = chunk[3]
                
                # Match the read ID
                match = fqIDre.match(header_line)
                if match:
                    currentID = match.group(1)  # Extract ID without '@'
                    if currentID in trimdict:
                        trim_amount = trimdict[currentID]
                        stats["trim_amounts"][trim_amount] += 1
                        
                        # Annotate the header line
                        if trim_amount > 0:
                            new_header = header_line.strip() + f"  {trim_amount} removed from left\n"
                            stats["trimmed_left"] += 1
                        elif trim_amount < 0:
                            new_header = header_line.strip() + f"  {abs(trim_amount)} removed from right\n"
                            stats["trimmed_right"] += 1
                        else:
                            new_header = header_line
                            stats["not_trimmed"] += 1
                        
                        # Process sequence line
                        if trim_amount > 0:
                            new_seq = seq_line[trim_amount:]
                            new_qual = qual_line[trim_amount:]
                        elif trim_amount < 0:
                            # Trim from the right
                            new_seq = seq_line[:trim_amount] + "\n"
                            new_qual = qual_line[:trim_amount] + "\n"
                        else:
                            new_seq = seq_line
                            new_qual = qual_line
                        
                        # Write the modified read
                        fastq_out.write(new_header)
                        fastq_out.write(new_seq)
                        fastq_out.write(plus_line)
                        fastq_out.write(new_qual)
                    else:
                        # No trim data found for this read
                        stats["not_trimmed"] += 1
                        fastq_out.write(header_line)
                        fastq_out.write(seq_line)
                        fastq_out.write(plus_line)
                        fastq_out.write(qual_line)
                else:
                    # No matching ID pattern
                    stats["no_match_id"] += 1
                    fastq_out.write(header_line)
                    fastq_out.write(seq_line)
                    fastq_out.write(plus_line)
                    fastq_out.write(qual_line)
                
                # Reset for next chunk
                chunk = []
    except Exception as e:
        log(f"ERROR during FASTQ processing: {str(e)}")
    finally:
        fastq_in.close()
        fastq_out.close()

    # Report statistics
    elapsed = time.time() - trim_start_time
    log(f"Trimming complete in {elapsed:.1f} seconds")
    log(f"Total reads processed: {stats['total_reads']:,}")
    log(f"Reads trimmed from left: {stats['trimmed_left']:,}")
    log(f"Reads trimmed from right: {stats['trimmed_right']:,}")
    log(f"Reads not trimmed: {stats['not_trimmed']:,}")
    log(f"Reads with no matching ID pattern: {stats['no_match_id']:,}")
    
    # Write detailed report if requested
    if args.report:
        try:
            with open(args.report, "w") as report_file:
                report_file.write("# Trimming Report\n")
                report_file.write(f"Input file: {args.fastq}\n")
                report_file.write(f"Output file: {args.out}\n")
                report_file.write(f"Trim data: {args.trimdata}\n")
                report_file.write(f"Date: {datetime.now()}\n\n")
                
                report_file.write("## Statistics\n")
                report_file.write(f"Total reads processed: {stats['total_reads']:,}\n")
                report_file.write(f"Reads trimmed from left: {stats['trimmed_left']:,}\n")
                report_file.write(f"Reads trimmed from right: {stats['trimmed_right']:,}\n")
                report_file.write(f"Reads not trimmed: {stats['not_trimmed']:,}\n")
                report_file.write(f"Reads with no matching ID pattern: {stats['no_match_id']:,}\n\n")
                
                report_file.write("## Most Common Trim Amounts\n")
                for amount, count in stats["trim_amounts"].most_common(20):
                    report_file.write(f"{amount}: {count:,}\n")
            
            log(f"Detailed report written to {args.report}")
        except IOError as e:
            log(f"ERROR writing report: {str(e)}")

if __name__ == "__main__":
    main()
