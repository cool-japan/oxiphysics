#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

echo "Building oxiphysics WASM package..."
wasm-pack build --target web --out-name oxiphysics --scope cooljapan
echo "Done. Output in pkg/"
