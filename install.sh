#!/usr/bin/env bash
set -euo pipefail

# -------------------------
# ME-ICA source / runtime release selection
# -------------------------
# The source release/tag and runtime release/tag can be set independently.
# This is useful for staging or interim source releases that should use an
# already-published compatible runtime, for example:
#
#   MEICA_SOURCE_TAG="4.0.1-stage" \
#   MEICA_RUNTIME_TAG="4.0.1" \
#   bash install.sh "$HOME"
#
# Backward compatibility: MEICA_VERSION is treated as an alias for
# MEICA_SOURCE_TAG if MEICA_SOURCE_TAG is not set.
MEICA_SOURCE_TAG="${MEICA_SOURCE_TAG:-${MEICA_VERSION:-4.0.1}}"
MEICA_RUNTIME_REPO="${MEICA_RUNTIME_REPO:-https://github.com/ME-ICA/meica-afni}"
MEICA_RUNTIME_TAG="${MEICA_RUNTIME_TAG:-${MEICA_RUNTIME_VERSION:-${MEICA_SOURCE_TAG}}}"
MEICA_RUNTIME_ASSET="${MEICA_RUNTIME_ASSET:-meica-afni-runtime-linux-x86_64.tar.gz}"

# By default, the source checkout directory and micromamba environment name
# follow the source tag. Override these only if you need a custom local name.
MEICA_INSTALL_LABEL="${MEICA_INSTALL_LABEL:-${MEICA_SOURCE_TAG}}"
MEICA_ENV_NAME="${MEICA_ENV_NAME:-meicapy${MEICA_INSTALL_LABEL}}"

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
STARTDIR="$(pwd)"
if [ "${INSTALLDIR:0:1}" != "/" ]; then
  INSTALLDIR="$STARTDIR/$INSTALLDIR"
fi

MEICA_SRC_DIR="$INSTALLDIR/me-ica-${MEICA_INSTALL_LABEL}"

# -------------------------
# Platform detection
# -------------------------
OS_NAME="$(uname -s)"
ARCH_NAME="$(uname -m)"
IS_LINUX=0
IS_MAC=0
IS_WSL=0
MICROMAMBA_PLATFORM=""

case "$OS_NAME" in
  Linux)
    IS_LINUX=1
    MICROMAMBA_PLATFORM="linux-64"
    if grep -qi microsoft /proc/version 2>/dev/null; then
      IS_WSL=1
    fi
    ;;
  Darwin)
    IS_MAC=1
    case "$ARCH_NAME" in
      arm64|aarch64)
        MICROMAMBA_PLATFORM="osx-arm64"
        ;;
      x86_64)
        MICROMAMBA_PLATFORM="osx-64"
        ;;
      *)
        echo "Unsupported macOS architecture: $ARCH_NAME"
        exit 1
        ;;
    esac
    ;;
  *)
    echo "Unsupported operating system: $OS_NAME"
    echo "ME-ICA installation currently supports Linux/WSL2 and macOS."
    exit 1
    ;;
esac

# -------------------------
# macOS AFNI policy
# -------------------------
# A distributed macOS AFNI runtime is not yet provided by this installer.
# If installing on macOS, require an existing AFNI installation on PATH.
if [ "$IS_MAC" -eq 1 ]; then
  if ! command -v 3dinfo >/dev/null 2>&1 || ! command -v 3dSkullStrip >/dev/null 2>&1; then
    echo
    echo "macOS detected, but no distributed macOS ME-ICA AFNI runtime is available yet."
    echo
    echo "Please install AFNI for macOS and ensure AFNI commands are on PATH before"
    echo "running this installer. At minimum, this installer expects commands such as:"
    echo "  3dinfo"
    echo "  3dSkullStrip"
    echo
    echo "After AFNI is available on PATH, re-run:"
    echo "  bash install.sh $INSTALLDIR"
    echo
    exit 1
  fi
fi

# -------------------------
# Check for existing environment
# -------------------------
ENV_ROOT="$INSTALLDIR"
ENV_NAME="${MEICA_ENV_NAME}"
ENV_DIR="$ENV_ROOT/envs/$ENV_NAME"
if [ -e "$ENV_DIR" ]; then
  echo "Environment already exists at $ENV_DIR . Remove and re-run installer."
  exit 1
fi

# -------------------------
# Clone repo
# -------------------------
mkdir -p "$INSTALLDIR"
cd "$INSTALLDIR"
if [ -e "$MEICA_SRC_DIR" ]; then
  echo "ME-ICA source directory already exists at $MEICA_SRC_DIR"
  echo "Remove it and re-run installer."
  exit 1
fi

git clone -b "${MEICA_SOURCE_TAG}" --single-branch https://github.com/ME-ICA/me-ica.git "$MEICA_SRC_DIR"

# -------------------------
# Download distributed AFNI runtime on Linux/WSL
# -------------------------
# The runtime is installed inside the matched ME-ICA source tree:
#
#   $INSTALLDIR/
#     me-ica-${MEICA_INSTALL_LABEL}/
#       meica.py
#       libmeica/
#       meica-afni-runtime/
#
# This matches libmeica/runtime.py, which searches first MEICA_AFNI_RUNTIME and
# then the sibling runtime directory inside the ME-ICA source tree.
RUNTIME_INSTALLED=0
MEICA_AFNI_RUNTIME="$MEICA_SRC_DIR/meica-afni-runtime"
MEICA_AFNI_RUNTIME_ACTIVATE="$MEICA_AFNI_RUNTIME/activate_meica_afni"
MEICA_RUNTIME_UTILS="$MEICA_AFNI_RUNTIME/utils"

if [ "$IS_LINUX" -eq 1 ]; then
  echo
  echo "Downloading ME-ICA distributed AFNI runtime release:"
  echo "  source tag:  $MEICA_SOURCE_TAG"
  echo "  source dir:  $MEICA_SRC_DIR"
  echo "  runtime repo: $MEICA_RUNTIME_REPO"
  echo "  runtime tag:  $MEICA_RUNTIME_TAG"
  echo "  runtime file: $MEICA_RUNTIME_ASSET"
  echo

  RUNTIME_URL="${MEICA_RUNTIME_REPO}/releases/download/${MEICA_RUNTIME_TAG}/${MEICA_RUNTIME_ASSET}"
  RUNTIME_TARBALL="$INSTALLDIR/${MEICA_RUNTIME_ASSET}"

  if [ -e "$MEICA_AFNI_RUNTIME" ]; then
    echo "Runtime directory already exists at $MEICA_AFNI_RUNTIME"
    echo "Remove it and re-run installer."
    exit 1
  fi

  curl -L --fail --retry 3 --retry-delay 2 \
    -o "$RUNTIME_TARBALL" \
    "$RUNTIME_URL"

  # Extract into the matched ME-ICA source directory, so the runtime appears at:
  #   me-ica-${MEICA_INSTALL_LABEL}/meica-afni-runtime
  tar -xzf "$RUNTIME_TARBALL" -C "$MEICA_SRC_DIR"

  if [ ! -d "$MEICA_AFNI_RUNTIME" ]; then
    echo "ERROR: Runtime tarball did not create expected directory:"
    echo "  $MEICA_AFNI_RUNTIME"
    echo
    echo "Expected the tarball to contain a top-level directory named:"
    echo "  meica-afni-runtime"
    exit 1
  fi

  if [ ! -f "$MEICA_AFNI_RUNTIME_ACTIVATE" ]; then
    echo "ERROR: Missing runtime activation script:"
    echo "  $MEICA_AFNI_RUNTIME_ACTIVATE"
    exit 1
  fi

  if [ ! -x "$MEICA_RUNTIME_UTILS/dcm2niix" ]; then
    echo "ERROR: Runtime was installed, but dcm2niix was not found at:"
    echo "  $MEICA_RUNTIME_UTILS/dcm2niix"
    exit 1
  fi

  RUNTIME_INSTALLED=1
fi

# -------------------------
# Install micromamba locally
# -------------------------
mkdir -p "$INSTALLDIR/micromamba"
curl -Ls "https://micro.mamba.pm/api/micromamba/${MICROMAMBA_PLATFORM}/latest" \
  | tar -xvj -C "$INSTALLDIR/micromamba" bin/micromamba

export PATH="$INSTALLDIR/micromamba/bin:$PATH"
export MAMBA_ROOT_PREFIX="$ENV_ROOT"

# -------------------------
# Initialize micromamba shell
# -------------------------
eval "$("$INSTALLDIR/micromamba/bin/micromamba" shell hook --shell bash)"

# -------------------------
# Create environment
# -------------------------
# Do not install conda AFNI on Linux/WSL when the distributed runtime is present.
# On macOS, AFNI is expected to be provided by the user's system PATH.
micromamba create -y -n "$ENV_NAME" -c conda-forge \
  --strict-channel-priority \
  python=3.12.3 pip tcsh

micromamba activate "$ENV_NAME"

# Make runtime variables visible for the installer-time dependency check.
if [ "$RUNTIME_INSTALLED" -eq 1 ]; then
  export MEICA_AFNI_RUNTIME="$MEICA_AFNI_RUNTIME"
  export MEICA_AFNI_RUNTIME_ACTIVATE="$MEICA_AFNI_RUNTIME_ACTIVATE"
  export MEICA_REQUIRE_DISTRIBUTED_AFNI=1
  export PATH="$MEICA_RUNTIME_UTILS:$PATH"
else
  export MEICA_REQUIRE_DISTRIBUTED_AFNI=0
fi

export PATH="$MEICA_SRC_DIR:$PATH"

# -------------------------
# Install dependencies
# -------------------------
pip install -r "$MEICA_SRC_DIR/requirements.txt"

# -------------------------
# Create activation helper
# -------------------------
if [ "$RUNTIME_INSTALLED" -eq 1 ]; then
cat > "$HOME/activate_meica" <<EOF
# Source this file to activate ME-ICA
export PATH="$INSTALLDIR/micromamba/bin:\$PATH"
eval "\$(micromamba shell hook --shell bash)"
export MAMBA_ROOT_PREFIX=$ENV_ROOT
micromamba activate $ENV_NAME

# ME-ICA source
export MEICA_REPO="$MEICA_SRC_DIR"
export MEICA_SOURCE_TAG="$MEICA_SOURCE_TAG"
export MEICA_RUNTIME_TAG="$MEICA_RUNTIME_TAG"
export PATH="\$MEICA_REPO:\$PATH"

# ME-ICA distributed AFNI runtime.
# The runtime is installed inside the matched ME-ICA source tree:
#   \$MEICA_REPO/meica-afni-runtime
#
# The full runtime/bin directory is not placed on PATH globally.
# meica.py activates the runtime internally for dependency checks and processing.
export MEICA_AFNI_RUNTIME="\$MEICA_REPO/meica-afni-runtime"
export MEICA_AFNI_RUNTIME_ACTIVATE="\$MEICA_AFNI_RUNTIME/activate_meica_afni"
export MEICA_REQUIRE_DISTRIBUTED_AFNI=1

# Curated user-facing utilities from the runtime.
# This exposes dcm2niix without exposing the full AFNI runtime globally.
export PATH="\$MEICA_AFNI_RUNTIME/utils:\$PATH"

deactivate() {
#  while [ "\${CONDA_SHLVL:-0}" -gt 0 ]; do
    micromamba deactivate
#  done
}

meica.py

echo "Welcome to the ME-ICA environment"
echo
echo "ME-ICA is active for this shell."
echo
echo "For help:"
echo "  meica.py -h"
echo
echo "To analyze data from the data directory:"
echo "  cd to_datadir"
echo "  meica.py -d data1_e*nii.gz"
echo
echo "To convert a folder of DICOMs to NIfTI with dcm2niix:"
echo "  mkdir -p nifti_out"
echo "  dcm2niix -z y -b y -o nifti_out /path/to/dicom_folder"
echo
echo "To deactivate environment:"
echo "  deactivate"
EOF
else
cat > "$HOME/activate_meica" <<EOF
# Source this file to activate ME-ICA
export PATH="$INSTALLDIR/micromamba/bin:\$PATH"
eval "\$(micromamba shell hook --shell bash)"
export MAMBA_ROOT_PREFIX=$ENV_ROOT
micromamba activate $ENV_NAME

# ME-ICA source
export MEICA_REPO="$MEICA_SRC_DIR"
export MEICA_SOURCE_TAG="$MEICA_SOURCE_TAG"
export MEICA_RUNTIME_TAG="$MEICA_RUNTIME_TAG"
export PATH="\$MEICA_REPO:\$PATH"

# No distributed AFNI runtime is installed on this platform by this installer.
# AFNI is expected to be available from the user's system PATH.
export MEICA_REQUIRE_DISTRIBUTED_AFNI=0

deactivate() {
#  while [ "\${CONDA_SHLVL:-0}" -gt 0 ]; do
    micromamba deactivate
#  done
}

meica.py

echo "Welcome to the ME-ICA environment"
echo
echo "ME-ICA is active for this shell."
echo
echo "For help:"
echo "  meica.py -h"
echo
echo "To analyze data from the data directory:"
echo "  cd to_datadir"
echo "  meica.py -d data1_e*nii.gz"
echo
echo "This platform is using AFNI from the system PATH."
echo
echo "To deactivate environment:"
echo "  deactivate"
EOF
fi

# -------------------------
# Test install
# -------------------------
cd "$STARTDIR"
set +e
"$MEICA_SRC_DIR/meica.py"
set -e

# -------------------------
# Installation message
# -------------------------
echo
echo "Installation complete. Review dependency check from meica.py (above)."
echo
echo "Installed source tag:"
echo "  $MEICA_SOURCE_TAG"
echo "Installed source directory:"
echo "  $MEICA_SRC_DIR"
if [ "$RUNTIME_INSTALLED" -eq 1 ]; then
  echo "Installed runtime tag:"
  echo "  $MEICA_RUNTIME_TAG"
fi
echo
echo "To begin using ME-ICA:"
echo "  source ~/activate_meica"
echo
echo "For help:"
echo "  meica.py -h"
echo
echo "To analyze data:"
echo "  cd to_datadir"
echo "  meica.py -d data1_e*nii.gz"
if [ "$RUNTIME_INSTALLED" -eq 1 ]; then
  echo
  echo "This installation includes the distributed ME-ICA AFNI runtime:"
  echo "  $MEICA_AFNI_RUNTIME"
  echo
  echo "The runtime is used internally by meica.py and generated _meica_*.sh scripts."
  echo "The full AFNI runtime is not globally exposed on PATH."
  echo
  echo "Selected runtime utilities are exposed when the environment is activated."
  echo "For DICOM conversion after running 'source ~/activate_meica':"
  echo "  mkdir -p nifti_out"
  echo "  dcm2niix -z y -b y -o nifti_out /path/to/dicom_folder"
else
  echo
  echo "This platform is using AFNI from the system PATH."
fi
echo
echo "To deactivate environment:"
echo "  deactivate"

deactivate() {
#  while [ "${CONDA_SHLVL:-0}" -gt 0 ]; do
    micromamba deactivate
#  done
}
deactivate
