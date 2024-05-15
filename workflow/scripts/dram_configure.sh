#!/usr/bin/env bash

set -euo pipefail

db_dir="/home/jlanga/share/db/dram"

DRAM-setup.py prepare_databases \
    --output_dir $db_dir
