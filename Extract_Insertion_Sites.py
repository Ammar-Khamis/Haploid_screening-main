#!/usr/bin/env python3
"""
Extract_Insertion_Sites.py
-------------------------
Extracts insertion sites from SAM/BAM file by looking at the 5' end
of aligned reads, accounting for the strand.

For PiggyBac transposon insertions, the 5' end of the read should 
correspond to the TTAA insertion site.

Input: SAM format from stdin (e.g., output of `samtools view`)
Output: BED format with insertion sites

Usage Example:
  samtools view mapped.bam | python Extract_Insertion_Sites.py --output insertion_sites.bed
"""
import sys
import argparse
import re
import time
from datetime import datetime
from collections import Counter, defaultdict

def log(message):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"[{timestamp}] {message}", file=sys.stderr)

def parse_cigar(cigar_str):
    """Parse CIGAR string and return operations"""
    cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar_str)
    return [(int(count), op) for count, op in cigar_ops]

def calculate_5prime_position(pos, cigar, is_reverse):
    """
    Calculate the 5' position based on alignment position and CIGAR string
    
    Args:
        pos: 1-based leftmost position from SAM
        cigar: Parsed CIGAR operations
        is_reverse: True if read aligned to reverse strand
        
    Returns:
        The 5' genomic position (0-based for BED output)
    """
    pos_0based = pos - 1  # Convert to 0-based for BED
    
    if is_reverse:
        # For reverse strand, 5' is at the end of the alignment
        # Need to walk through CIGAR to find the rightmost position
        for count, op in cigar:
            if op in 'MDN=X':  # These operations consume reference
                pos_0based += count
        
        # Adjust to get the actual 5' position (subtract 1 for 0-based end position)
        return pos_0based - 1
    else:
        # For forward strand, 5' is at the start of the alignment
        # No adjustment needed for 0-based BED start
        return pos_0based

def main():
    parser = argparse.ArgumentParser(
        description="Extract insertion sites from SAM/BAM data."
    )
    parser.add_argument("--output", required=True,
        help="Output BED file for insertion sites.")
    parser.add_argument("--min_mapq", type=int, default=10,
        help="Minimum mapping quality to include (default: 10)")
    parser.add_argument("--merge_window", type=int, default=5,
        help="Window size to merge nearby insertions (default: 5bp)")

    args = parser.parse_args()
    start_time = time.time()
    
    # Open output file
    try:
        out_file = open(args.output, "w")
    except IOError as e:
        log(f"ERROR opening output file: {str(e)}")
        sys.exit(1)

    log(f"Extracting insertion sites to {args.output}")
    log(f"Minimum MAPQ: {args.min_mapq}")
    
    # Statistics tracking
    stats = {
        "total_reads": 0,
        "passed_mapq": 0,
        "failed_mapq": 0,
        "forward_strand": 0,
        "reverse_strand": 0,
        "chromosomes": Counter()
    }
    
    # Store insertion sites to merge nearby ones
    insertion_sites = defaultdict(list)  # chrom -> [(pos, read_id), ...]
    
    # Process SAM input from stdin
    try:
        for line_num, line in enumerate(sys.stdin):
            # Skip header lines
            if line.startswith('@'):
                continue
            
            stats["total_reads"] += 1
            
            # Progress indicator
            if stats["total_reads"] % 100000 == 0:
                elapsed = time.time() - start_time
                log(f"Processed {stats['total_reads']:,} alignments ({elapsed:.1f} seconds)")
            
            # Parse SAM fields
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
                
            read_id = fields[0]
            flag = int(fields[1])
            chrom = fields[2]
            pos = int(fields[3])  # 1-based leftmost position
            mapq = int(fields[4])
            cigar_str = fields[5]
            
            # Skip unmapped reads
            if (flag & 0x4) != 0 or chrom == '*':
                continue
                
            # Apply MAPQ filter
            if mapq < args.min_mapq:
                stats["failed_mapq"] += 1
                continue
                
            stats["passed_mapq"] += 1
            stats["chromosomes"][chrom] += 1
            
            # Determine strand
            is_reverse = (flag & 0x10) != 0
            if is_reverse:
                stats["reverse_strand"] += 1
                strand = '-'
            else:
                stats["forward_strand"] += 1
                strand = '+'
            
            # Parse CIGAR and calculate 5' position
            cigar_ops = parse_cigar(cigar_str)
            site_pos = calculate_5prime_position(pos, cigar_ops, is_reverse)
            
            # Store the insertion site
            insertion_sites[chrom].append((site_pos, read_id, strand))
    
    except Exception as e:
        log(f"ERROR processing SAM input: {str(e)}")
    
    # Merge nearby insertion sites and write to BED file
    log("Merging nearby insertion sites...")
    total_sites = 0
    merged_sites = 0
    
    for chrom, sites in insertion_sites.items():
        # Sort sites by position
        sites.sort()
        
        # Process sites to merge those within the window
        current_start = -1
        current_end = -1
        current_ids = set()
        current_strands = []
        
        for pos, read_id, strand in sites:
            total_sites += 1
            
            if current_start == -1:
                # First site
                current_start = pos
                current_end = pos + 1  # BED end is exclusive
                current_ids.add(read_id)
                current_strands.append(strand)
            elif pos <= current_start + args.merge_window:
                # Merge with current site
                current_end = pos + 1  # Update end position
                current_ids.add(read_id)
                current_strands.append(strand)
                merged_sites += 1
            else:
                # Write current site and start a new one
                count = len(current_ids)
                # Determine dominant strand
                strand_counts = Counter(current_strands)
                dominant_strand = '+' if strand_counts['+'] >= strand_counts['-'] else '-'
                
                out_file.write(f"{chrom}\t{current_start}\t{current_end}\t{count}\t{count}\t{dominant_strand}\n")
                
                # Start new site
                current_start = pos
                current_end = pos + 1
                current_ids = {read_id}
                current_strands = [strand]
        
        # Write the last site if there is one
        if current_start != -1:
            count = len(current_ids)
            # Determine dominant strand
            strand_counts = Counter(current_strands)
            dominant_strand = '+' if strand_counts['+'] >= strand_counts['-'] else '-'
            
            out_file.write(f"{chrom}\t{current_start}\t{current_end}\t{count}\t{count}\t{dominant_strand}\n")
    
    out_file.close()
    
    # Report statistics
    elapsed = time.time() - start_time
    log(f"Extraction complete in {elapsed:.1f} seconds")
    log(f"Total reads processed: {stats['total_reads']:,}")
    log(f"Reads passing MAPQ filter: {stats['passed_mapq']:,} ({(stats['passed_mapq']/stats['total_reads']*100 if stats['total_reads'] > 0 else 0):.1f}%)")
    log(f"Forward strand reads: {stats['forward_strand']:,}")
    log(f"Reverse strand reads: {stats['reverse_strand']:,}")
    log(f"Total insertion sites: {total_sites:,}")
    log(f"Merged sites: {merged_sites:,}")
    log(f"Final unique sites: {total_sites - merged_sites:,}")
    
    # Top chromosomes
    log("Top chromosomes with insertions:")
    for chrom, count in stats["chromosomes"].most_common(10):
        log(f"  {chrom}: {count:,}")

if __name__ == "__main__":
    main()
