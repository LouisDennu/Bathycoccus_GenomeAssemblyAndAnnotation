#!/bin/bash
set -euo pipefail

TEMPLATE="singularity/singu.def"
ENV_DIR="singularity/env_yml"
OUT_DIR="singularity"
mkdir -p "$OUT_DIR"

for yml in "$ENV_DIR"/*.yml; do
    base=$(basename "$yml" .yml)       # e.g. flye, medaka, pilon...
    env_name=$(grep '^name:' "$yml" | awk '{print $2}')
    def_file="$OUT_DIR/${base}.def"
    sif_file="$OUT_DIR/${base}.sif"

    echo "Generating $def_file and building $sif_file..."

    # Replace placeholders in template
    sed "s|ENV_YML|$yml|g; s|ENV_NAME|$env_name|g" "$TEMPLATE" > "$def_file"

    # Build the image with desired name
    singularity build --fakeroot "$sif_file" "$def_file"
done
