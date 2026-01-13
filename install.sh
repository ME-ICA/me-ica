#!/usr/bin/env bash
set -euo pipefail
MEICA_VERSION=3.9.7

# -------------------------
# Check usage
# -------------------------
INSTALLDIR="${1:-}"
if [ -z "$INSTALLDIR" ]; then
  echo "Usage: install.sh INSTALLATION_DIRECTORY"
  exit 1
fi

# -------------------------
# Setup directories
# -------------------------
mkdir -p "$INSTALLDIR"
cd "$INSTALLDIR"

# -------------------------
# Clone repo
# -------------------------
git clone -b ${MEICA_VERSION}-gpl --single-branch https://github.com/ME-ICA/me-ica.git me-ica-${MEICA_VERSION}

# -------------------------
# Install micromamba locally
# -------------------------
mkdir -p "$INSTALLDIR/micromamba"
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C "$INSTALLDIR/micromamba" bin/micromamba

export PATH="$INSTALLDIR/micromamba/bin:$PATH"

# -------------------------
# Initialize micromamba shell
# -------------------------
eval "$("$INSTALLDIR/micromamba/bin/micromamba" shell hook --shell bash)"

# -------------------------
# Create environment
# -------------------------
micromamba create -y -n meicapy -c conda-forge \
  --strict-channel-priority \
  python=3.12 pip afni gsl
micromamba activate meicapy

# -------------------------
# Install dependencies
# -------------------------
pip install -r "$INSTALLDIR/me-ica-${MEICA_VERSION}/requirements.txt"

# -------------------------
# Test install
# -------------------------
"$INSTALLDIR/me-ica-${MEICA_VERSION}/meica.py"

# -------------------------
# Create activation helper
# -------------------------
cat > "$HOME/activate_meica" <<EOF
# Source this file to activate ME-ICA
export PATH="$INSTALLDIR/micromamba/bin:\$PATH"
eval "\$(micromamba shell hook --shell bash)"
micromamba activate meicapy
export PATH="$INSTALLDIR/me-ica-${MEICA_VERSION}/:\$PATH"
deactivate() {
  while [ "${CONDA_SHLVL:-0}" -gt 0 ]; do
    micromamba deactivate
  done
}
echo "Welcome to the ME-ICA environment"
echo "To use (from the data directory!):"
echo "  meica.py -h"
echo "To exit:"
echo "  deactivate"
EOF

echo
echo "Installation complete."
echo "To begin using ME-ICA:"
echo "  source ~/activate_meica"
echo "For help:"
echo "  meica.py -h"
echo "To analyze data:"
echo "  cd to_datadir"
echo "  meica.py -d data1_e*nii.gz"
echo "To deactivate environment:"
echo "  deactivate"

deactivate() {
  while [ "${CONDA_SHLVL:-0}" -gt 0 ]; do
    micromamba deactivate
  done
}
deactivate
