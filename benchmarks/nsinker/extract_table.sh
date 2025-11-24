#!/usr/bin/env bash
#
# extract_table.sh
#   Extract the 4‑column Section‑timings table from a text file
#   and convert it to CSV (Section,no. calls,wall time,% of total).
#
# Usage: ./extract_table.sh <input.txt> <output.csv>

if [ $# -ne 2 ]; then
  echo "Usage: $0 <input.txt> <output.csv>"
  exit 1
fi

INPUT="$1"
OUTPUT="$2"

# 1) sed -n '/start/,/end/p'  → grab from first 4‑col border to next one
# 2) sed '1d;$d'            → delete the two border lines we no longer need
# 3) sed 's/^| *//; s/ *| *$//; s/ *| */,/g'
#    • strip leading "| "
#    • strip trailing " |"
#    • turn " | " between fields into commas
#
sed -n '/^\+[-]*\+[-]*\+[-]*\+[-]*\+$/,/^\+[-]*\+[-]*\+[-]*\+[-]*\+$/p' "$INPUT" \
  | sed '1d;$d' \
  | sed 's/^| *//; s/ *| *$//; s/ *| */,/g' \
  > "$OUTPUT"

echo "✔ Parsed table written to $OUTPUT"
