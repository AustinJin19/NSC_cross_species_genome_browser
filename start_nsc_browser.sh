#!/usr/bin/env bash
set -euo pipefail

cd /data/upload/CS_Genome_Browser/NSC_Genome_Browser
source .venv/bin/activate
exec python cross_species_visual_tool/web_app/server.py --port 8765
