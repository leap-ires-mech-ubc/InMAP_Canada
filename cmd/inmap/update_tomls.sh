#!/usr/bin/env bash
# Updates OutputFile date/suffix, normalizes EmissionsShapefiles basenames,
# and checks file existence for configGEMMACH_scenario{1..9}.toml
# Usage: bash update_tomls.sh (Made with help of 365 Copilot)
set -u  # (avoid -e to allow collecting all missing files across inputs)
shopt -s nullglob

TODAY=20260223 #"$(date +%Y%m%d)"                     # For OutputFile
EMIS_DATE="${EMIS_DATE:-202600223}"          # For EmissionsShapefiles prefix
TOML_DIR="/home/tfmrodge/projects/def-rscholes/tfmrodge/InMAP_Canada/cmd/inmap/"

# New: optional base directory to apply to EmissionsShapefiles
# If empty, keep original directories. If set, override directory for rewritten items.
EMIS_DIR="/home/tfmrodge/scratch/GEMMACH_data/data/emissions_nc/20260223/"

# New: whether to also move 'majorpts' files into EMIS_DIR (keep their names)
# 0 = leave majorpts completely unchanged (default) #"${APPLY_DIR_TO_MAJORPTS:-0}"
# 1 = change only their directory to EMIS_DIR (basename unchanged)
APPLY_DIR_TO_MAJORPTS=0 

status=0

echo "==> Using TODAY=${TODAY} for OutputFile, EMIS_DATE=${EMIS_DATE} for EmissionsShapefiles"
if [[ -n "$EMIS_DIR" ]]; then
  echo "==> Overriding EmissionsShapefiles directory with: $EMIS_DIR"
  echo "==> APPLY_DIR_TO_MAJORPTS=$APPLY_DIR_TO_MAJORPTS (0=leave majorpts as-is, 1=move into EMIS_DIR)"
fi

for f in "${TOML_DIR}"configGEMMACH_scenario{1..9}.toml; do
  [[ -f "$f" ]] || continue

  echo
  echo "Processing: $f"

  tmp="$(mktemp)"
  emis_log="$(mktemp)"

  awk -v TODAY="$TODAY" \
      -v EMIS_DATE="$EMIS_DATE" \
      -v EMIS_LOG="$emis_log" \
      -v EMIS_DIR="$EMIS_DIR" \
      -v APPLY_DIR_TO_MAJORPTS="$APPLY_DIR_TO_MAJORPTS" '
    function ensure_trailing_slash(s) {
      if (s == "") return s
      return (s ~ /\/$/) ? s : s "/"
    }

    function rewrite_outputfile_path(p,   d,b,noext,new) {
      d = p; sub(/[^/]*$/,"",d)         # dir (keep trailing / if any)
      b = p; sub(/^.*\//,"",b)          # base
      noext = b; sub(/\.shp$/,"",noext) # strip .shp
      sub(/_static$/,"",noext)          # strip trailing _static (if present)
      sub(/^[0-9]{8}_/,"",noext)        # strip leading YYYYMMDD_ (if present)
      new = d TODAY "_" noext "_static.shp"
      return new
    }

    function rewrite_emis_path(p,   d,b,noext,is_majorpts,dir_new,new_b,newp) {
      d = p; sub(/[^/]*$/,"",d)                   # original dir
      b = p; sub(/^.*\//,"",b)                    # base (filename only)
      is_majorpts = (tolower(b) ~ /majorpts/)

      # Determine the directory to use
      dir_new = d
      if (EMIS_DIR != "") {
        dir_new = ensure_trailing_slash(EMIS_DIR)
      }

      # If "majorpts" and we are not applying new dir to it, leave completely unchanged
      if (is_majorpts && APPLY_DIR_TO_MAJORPTS == 0) {
        return p
      }

      # If "majorpts" and APPLY_DIR_TO_MAJORPTS==1: only rebase directory, keep filename
      if (is_majorpts && APPLY_DIR_TO_MAJORPTS == 1) {
        return dir_new b
      }

      # Non-majorpts: rewrite filename to EMIS_DATE_<core>_area_EPSG4326.shp
      noext = b; sub(/\.shp$/,"",noext)           # strip .shp
      sub(/^[0-9]{8}_/,"",noext)                  # strip leading YYYYMMDD_
      sub(/_[^_]+$/, "", noext)                   # drop last token (e.g., _majorpts/_areaEPSG4326/etc)
      new_b = EMIS_DATE "_" noext "_area_EPSG4326.shp"
      newp = dir_new new_b
      return newp
    }

    BEGIN { in_emis = 0 }

    {
      # Do not touch commented lines
      if ($0 ~ /^[[:space:]]*#/) { print; next }

      # --- OutputFile rewrite ---
      if ($0 ~ /^[[:space:]]*OutputFile[[:space:]]*=/) {
        line = $0
        if (match(line, /"[^"]*"/)) {
          quoted = substr(line, RSTART, RLENGTH)
          path = quoted; gsub(/^"/,"",path); gsub(/"$/,"",path)
          newp = rewrite_outputfile_path(path)
          $0 = substr(line,1,RSTART-1) "\"" newp "\"" substr(line, RSTART+RLENGTH)
        }
        print; next
      }

      # --- EmissionsShapefiles block detection ---
      if ($0 ~ /^[[:space:]]*EmissionsShapefiles[[:space:]]*=/) {
        in_emis = 1
      }

      if (in_emis) {
        # Replace any number of quoted paths in this line
        line = $0
        out = ""
        while (match(line, /"[^"]*"/)) {
          prefix = substr(line, 1, RSTART-1)
          quoted = substr(line, RSTART, RLENGTH)
          path = quoted; gsub(/^"/,"",path); gsub(/"$/,"",path)
          newp = rewrite_emis_path(path)
          out = out prefix "\"" newp "\""
          if (EMIS_LOG != "") { print newp >> EMIS_LOG }
          line = substr(line, RSTART + RLENGTH)
        }
        $0 = out line
        if ($0 ~ /\]/) { in_emis = 0 }
        print; next
      }

      # Default: print unchanged
      print
    }
  ' "$f" > "$tmp"

  # Backup and replace atomically
  cp -p "$f" "$f.bak"
  mv "$tmp" "$f"

  echo "  - Updated OutputFile date/suffix and normalized EmissionsShapefiles."

  # ---- Check existence of all referenced EmissionsShapefiles ----
  if [[ -s "$emis_log" ]]; then
    echo "  - Checking EmissionsShapefiles existence:"
    missing=0
    # Unique paths
    while IFS= read -r p; do
      if [[ -e "$p" ]]; then
        echo "      ✓ $p"
      else
        echo "      ✗ MISSING: $p"
        ((missing++))
      fi
    done < <(sort -u "$emis_log")

    if (( missing > 0 )); then
      echo "  ! ${missing} missing shapefile(s) in $f"
      status=1
    else
      echo "  - All shapefiles exist."
    fi
  else
    echo "  - No EmissionsShapefiles found in $f (nothing to check)."
  fi

  rm -f "$emis_log"
done

echo
if (( status == 0 )); then
  echo "Done. All files processed successfully."
else
  echo "Completed with missing shapefiles. See marks above."
fi

exit $status