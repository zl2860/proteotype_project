#!/usr/bin/env bash
set -euo pipefail

# Backward-compatible wrapper.
# The UKB .bim files shown by the user use rsIDs in column 2, so the maintained
# extraction path is rsID-based.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/extract_ukb_st10_core_by_rsid_on_server.sh" "$@"
