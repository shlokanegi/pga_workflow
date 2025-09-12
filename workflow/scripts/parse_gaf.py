from sys import argv,stderr
import zstandard as zstd
from contextlib import contextmanager

MAP_Q_ID = 11

@contextmanager
def open_gaf(filename):
    if filename.endswith(".zst"):
        with zstd.open(filename, 'rt') as f:
            yield f
    else:
        with open(filename) as f:
            yield f

# def decompress_zst(input_file, output_file):
#     with open(input_file, 'rb') as ifh:

def print_only_primary_alignments(file):
    star_lines = tag_lines = no_tag_lines = tag_lines_selected = no_tag_lines_selected = tot_lines = 0
    with open_gaf(file) as f:
        for line in f:
            tot_lines += 1
            sp_line = line.strip().split('\t')
            if sp_line[4] == '*':
                star_lines += 1
                continue
            elif sp_line[1].startswith("M"):
                tag_lines += 1
                if int(sp_line[MAP_Q_ID+2]) >= 0:
                    tag_lines_selected += 1
                    print(line, end="")
            else:
                no_tag_lines += 1
                if int(sp_line[MAP_Q_ID]) >= 0:
                    no_tag_lines_selected += 1
                    print(line, end="")
    print(f"TOT_LINES: {tot_lines}\tSTAR_LINES: {star_lines}\tML/MF_TAG_LINES: {tag_lines}\tSELECTED_ML/MF_TAG_LINES: {tag_lines_selected}\tNO_TAG_LINES: {no_tag_lines}\tSELECTED_NO_TAG_LINES: {no_tag_lines_selected}\tSUM: {star_lines+tag_lines+no_tag_lines}", file=stderr)
if __name__ == "__main__":
    print_only_primary_alignments(argv[1])
