#!/usr/bin/env python3
"""
extract_table.py

Extracts a fixed-format ASCII timing table from a text file
and writes it out as a CSV file.

Usage:
    python extract_table.py input.txt output.csv
"""
import re
import csv
import argparse
import sys

def extract_table(input_path, output_path):
    # Regex to match the ASCII table border lines
    border_re = re.compile(r"^\+[\-]+\+[\-]+\+[\-]+\+[\-]+\+$")

    # Read all lines from the input file
    with open(input_path, 'r') as f:
        lines = f.readlines()

    # Find the first and second border lines
    start_idx = None
    end_idx = None
    for i, line in enumerate(lines):
        if border_re.match(line):
            if start_idx is None:
                start_idx = i
            else:
                end_idx = i
                break

    if start_idx is None or end_idx is None:
        print("Error: Could not find the table borders in the file.", file=sys.stderr)
        sys.exit(1)

    # Extract table body (excluding the border lines)
    body_lines = lines[start_idx + 1:end_idx]

    rows = []
    for line in body_lines:
        line = line.rstrip()  # Strip trailing newline and spaces
        # Skip empty rows
        if not line.strip():
            continue
        # Must start and end with |
        if not (line.startswith('|') and line.endswith('|')):
            continue
        # Remove leading and trailing '|' and split on '|'
        parts = [cell.strip() for cell in line.strip('|').split('|')]
        # We're interested only in rows with exactly 4 columns
        if len(parts) == 4:
            rows.append(parts)

    if not rows:
        print("Error: No table rows found after parsing.", file=sys.stderr)
        sys.exit(1)

    # Write to CSV
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Header row
        writer.writerow(rows[0])
        # Data rows
        for row in rows[1:]:
            writer.writerow(row)

    print(f"✔ Extracted {len(rows) - 1} data rows to '{output_path}'")


def main():
    parser = argparse.ArgumentParser(description='Extract ASCII timing table to CSV')
    parser.add_argument('input_file', help='Path to input text file')
    parser.add_argument('output_file', help='Path to output CSV file')
    args = parser.parse_args()
    extract_table(args.input_file, args.output_file)


if __name__ == '__main__':
    main()
