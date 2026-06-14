#!/usr/bin/env bash
set -euo pipefail
MEICA_VERSION=4.0.0

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
STARTDIR=`pwd`
if [ "${INSTALLDIR:0:1}" != "/" ]; then
	INSTALLDIR=$STARTDIR/$INSTALLDIR
fi


# -------------------------
# Check for existing environment
# -------------------------
ENV_ROOT="$INSTALLDIR"
ENV_NAME="meicapy${MEICA_VERSION}"
ENV_DIR="$ENV_ROOT/envs/$ENV_NAME"
if [ -e "$ENV_DIR" ]; then
	echo "Environment already exists at $ENV_DIR . Remove and re-run installer."
	exit
fi

# -------------------------
# Clone repo
# -------------------------
mkdir -p "$INSTALLDIR"
cd "$INSTALLDIR"
git clone -b ${MEICA_VERSION}-gpl --single-branch https://github.com/ME-ICA/me-ica.git me-ica-${MEICA_VERSION}

# -------------------------
# Install micromamba locally
# -------------------------
mkdir -p $INSTALLDIR/micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C "$INSTALLDIR/micromamba" bin/micromamba
export PATH="$INSTALLDIR/micromamba/bin:$PATH"
export MAMBA_ROOT_PREFIX=$ENV_ROOT

# -------------------------
# Initialize micromamba shell
# -------------------------
eval "$("$INSTALLDIR/micromamba/bin/micromamba" shell hook --shell bash)"

# -------------------------
# Create environment
# -------------------------
micromamba create -y -n $ENV_NAME -c conda-forge \
  --strict-channel-priority \
  python=3.12.3 pip afni gsl
micromamba activate $ENV_NAME

# -------------------------
# Install dependencies
# -------------------------
pip install -r "$INSTALLDIR/me-ica-${MEICA_VERSION}/requirements.txt"

# -------------------------
# Create activation helper
# -------------------------
cat > "$HOME/activate_meica" <<EOF
# Source this file to activate ME-ICA
export PATH="$INSTALLDIR/micromamba/bin:\$PATH"
eval "\$(micromamba shell hook --shell bash)"
export MAMBA_ROOT_PREFIX=$ENV_ROOT
micromamba activate $ENV_NAME
export PATH="$INSTALLDIR/me-ica-${MEICA_VERSION}/:\$PATH"
deactivate() {
#  while [ "\${CONDA_SHLVL:-0}" -gt 0 ]; do
    micromamba deactivate
#  done
}
meica.py
echo "Welcome to the ME-ICA environment"
echo "To use (from the data directory!):"
echo "  meica.py -h"
echo "To exit:"
echo "  deactivate"
EOF

# -------------------------
# Test install
# -------------------------
cd $STARTDIR
set +e
"$INSTALLDIR/me-ica-${MEICA_VERSION}/meica.py"
set -e

# -------------------------
# Test install
# -------------------------
echo
echo "Installation complete. Review dependency check from meica.py (above)."
echo
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
#  while [ "${CONDA_SHLVL:-0}" -gt 0 ]; do
    micromamba deactivate
#  done
}
deactivate
