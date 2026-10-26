#!/usr/bin/env bash
# Download the open BBBC039 dataset (~80 MB) into ./data/BBBC039
# Source: Broad Bioimage Benchmark Collection, https://bbbc.broadinstitute.org/BBBC039 (CC0)
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/BBBC039
cd data/BBBC039
for f in images masks metadata; do
  if [ ! -d "$f" ]; then
    echo "Downloading $f.zip ..."
    curl -L --fail -o "$f.zip" "https://data.broadinstitute.org/bbbc/BBBC039/$f.zip"
    unzip -q -o "$f.zip" && rm "$f.zip"
  fi
done
rm -rf __MACOSX
echo "Done: $(ls images | wc -l) images, $(ls masks | wc -l) masks in data/BBBC039"
