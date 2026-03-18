#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="${1:-tracked}"
MAX_SIZE_MB="${MAX_GIT_FILE_SIZE_MB:-50}"
MAX_SIZE_BYTES=$((MAX_SIZE_MB * 1024 * 1024))

collect_files() {
    case "$MODE" in
        --staged)
            git -C "$ROOT_DIR" diff --cached --name-only --diff-filter=AM -z
            ;;
        tracked)
            git -C "$ROOT_DIR" ls-files -z
            ;;
        *)
            echo "Usage: $0 [tracked|--staged]" >&2
            exit 2
            ;;
    esac
}

report_large_files() {
    local found=0
    local file
    local size

    while IFS= read -r -d '' file; do
        [ -f "$ROOT_DIR/$file" ] || continue
        size=$(stat -c '%s' "$ROOT_DIR/$file")
        if [ "$size" -gt "$MAX_SIZE_BYTES" ]; then
            if [ "$found" -eq 0 ]; then
                echo "Tracked file size policy violation:" >&2
            fi
            found=1
            printf '  %s (%s bytes)\n' "$file" "$size" >&2
        fi
    done

    if [ "$found" -ne 0 ]; then
        cat >&2 <<EOF
Files above ${MAX_SIZE_MB} MB should stay out of git in this repository.
Use ignored folders like data/output/ or data/input/local/, or external storage.
EOF
        exit 1
    fi
}

collect_files | report_large_files
