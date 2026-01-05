#!/usr/bin/env bash
# Install an RDP-trained classifier into classifiers/<marker>/<name>/ and update the classifier override config.
# Requires: bash, curl (for URLs), unzip/tar for archives, python3 with PyYAML.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: scripts/install_classifier.sh (--url URL | --file PATH) [--marker MARKER] [--name NAME] [--override PATH] [--force]

Options:
  --url URL        Remote archive containing rRNAClassifier.properties
  --file PATH      Local archive file (zip/tar.gz/tar.bz2) containing rRNAClassifier.properties
  --marker MARKER  Marker name used for directory placement (default: COI)
  --name NAME      Friendly name for this classifier (default: derived from archive or properties folder)
  --override PATH  Override YAML to write/update (default: config/classifier_override.yaml)
  --force          Overwrite existing target directory if it exists
  -h, --help       Show this help
EOF
  exit 1
}

URL=""
FILE=""
MARKER="COI"
NAME=""
OVERRIDE="config/classifier_override.yaml"
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --url) URL="${2:-}"; shift 2 ;;
    --file) FILE="${2:-}"; shift 2 ;;
    --marker) MARKER="${2:-}"; shift 2 ;;
    --name) NAME="${2:-}"; shift 2 ;;
    --override) OVERRIDE="${2:-}"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

if [[ -z "$URL" && -z "$FILE" ]]; then
  echo "Error: --url or --file is required." >&2
  usage
fi

if [[ -n "$URL" && -n "$FILE" ]]; then
  echo "Error: choose only one of --url or --file." >&2
  usage
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CLASS_ROOT="$REPO_ROOT/classifiers"
TMPDIR="$(mktemp -d)"
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

ARCHIVE_PATH="$TMPDIR/archive"
if [[ -n "$URL" ]]; then
  echo "Downloading classifier from $URL ..."
  curl -L "$URL" -o "$ARCHIVE_PATH"
else
  ARCHIVE_PATH="$(realpath "$FILE")"
fi

extract_dir="$TMPDIR/extracted"
mkdir -p "$extract_dir"

case "$ARCHIVE_PATH" in
  *.zip) unzip -q "$ARCHIVE_PATH" -d "$extract_dir" ;;
  *.tar.gz|*.tgz) tar -xzf "$ARCHIVE_PATH" -C "$extract_dir" ;;
  *.tar.bz2|*.tbz2) tar -xjf "$ARCHIVE_PATH" -C "$extract_dir" ;;
  *) echo "Unsupported archive type: $ARCHIVE_PATH" >&2; exit 1 ;;
esac

props_path="$(find "$extract_dir" -name 'rRNAClassifier.properties' -print -quit || true)"
if [[ -z "$props_path" ]]; then
  echo "rRNAClassifier.properties not found in archive." >&2
  exit 1
fi

props_dir="$(cd "$(dirname "$props_path")" && pwd)"

if [[ -z "$NAME" ]]; then
  NAME="$(basename "$props_dir" | tr ' ' '_' )"
fi

TARGET_DIR="$CLASS_ROOT/$MARKER/$NAME"
if [[ -e "$TARGET_DIR" && $FORCE -ne 1 ]]; then
  echo "Target $TARGET_DIR already exists. Use --force to overwrite." >&2
  exit 1
fi

rm -rf "$TARGET_DIR"
mkdir -p "$(dirname "$TARGET_DIR")"
cp -a "$props_dir" "$TARGET_DIR"

REL_PATH="classifiers/$MARKER/$NAME/rRNAClassifier.properties"

python3 - <<'PY'
import sys
from pathlib import Path

try:
    import yaml  # type: ignore
except ImportError as exc:  # pragma: no cover
    sys.stderr.write("PyYAML is required (pip install pyyaml)\n")
    sys.exit(1)

override_path = Path(sys.argv[1]).resolve()
classifier_path = sys.argv[2]

def deep_merge(base, new):
    for k, v in new.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_merge(base[k], v)
        else:
            base[k] = v
    return base

existing = {}
if override_path.exists():
    existing = yaml.safe_load(override_path.read_text()) or {}

update = {"modules": {"classification": {"rdp": {"t": classifier_path}}}}

merged = deep_merge(existing, update)
override_path.parent.mkdir(parents=True, exist_ok=True)
override_path.write_text(yaml.safe_dump(merged, sort_keys=False))

print(f"Updated override config at {override_path} with rdp.t={classifier_path}")
PY \
  "$REPO_ROOT/$OVERRIDE" \
  "$REL_PATH"

echo "Installed classifier to $TARGET_DIR"
echo "Config override written to $OVERRIDE"
echo "Ensure Snakefile_ESV loads overrides (already supported); run Snakemake with your chosen workflow."
