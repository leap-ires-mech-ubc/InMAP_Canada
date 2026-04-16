#!/usr/bin/env bash
TOML_DIR="/home/tfmrodge/projects/def-agiang01/tfmrodge/InMAP_Canada/cmd/inmap/"
for f in "${TOML_DIR}"configGEMMACH_scenario{1..9}.toml.bak; do
  orig="${f%.bak}"
  if [[ -f "$f" ]]; then
    echo "Restoring $orig from $f"
    cp -p "$f" "$orig"
  fi
done