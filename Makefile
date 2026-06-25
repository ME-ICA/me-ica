# Convenience targets for building and using a ME-ICA Docker image.
#
# The Docker image intentionally uses the same install.sh path as a normal user
# install.  install.sh creates /root/activate_meica inside the image, and the
# container entrypoint sources that activation helper before running commands.
#
# The resulting image has:
#   /opt/meica/me-ica-$(MEICA_VERSION)/
#   /opt/meica/me-ica-$(MEICA_VERSION)/meica-afni-runtime/
#   /root/activate_meica
#
# AFNI is not installed as a system package.  The distributed runtime is used by
# meica.py/runtime.py, and runtime/utils is exposed by activate_meica for tools
# such as dcm2niix.

MEICA_VERSION ?= 4.0.1
IMAGE ?= meica:$(MEICA_VERSION)
INSTALLDIR ?= /opt/meica
DATA ?= $(PWD)
DOCKER ?= docker

.PHONY: build test run shell meica dcm2niix clean build_dev run_dev exec_dev

build:
	$(DOCKER) build \
		--build-arg MEICA_VERSION=$(MEICA_VERSION) \
		--build-arg INSTALLDIR=$(INSTALLDIR) \
		-t $(IMAGE) .

test:
	$(DOCKER) run --rm $(IMAGE) meica.py -h
	$(DOCKER) run --rm $(IMAGE) dcm2niix -h >/dev/null

# Run the image with the current directory mounted as /data.
run:
	$(DOCKER) run --rm -it \
		-v "$(DATA):/data" \
		-w /data \
		$(IMAGE)

shell:
	$(DOCKER) run --rm -it \
		-v "$(DATA):/data" \
		-w /data \
		$(IMAGE) bash

# Example:
#   make meica ARGS='-d sub-01_task-rest_e*nii.gz'
meica:
	$(DOCKER) run --rm -it \
		-v "$(DATA):/data" \
		-w /data \
		$(IMAGE) meica.py $(ARGS)

# Example:
#   make dcm2niix DICOM=/data/dicoms OUT=/data/nifti_out
DICOM ?= /data/dicom_folder
OUT ?= /data/nifti_out
dcm2niix:
	$(DOCKER) run --rm -it \
		-v "$(DATA):/data" \
		-w /data \
		$(IMAGE) bash -lc 'mkdir -p "$(OUT)" && dcm2niix -z y -b y -o "$(OUT)" "$(DICOM)"'

clean:
	-$(DOCKER) image rm $(IMAGE)

# Backward-compatible aliases for older development target names.
build_dev: build
run_dev: run
exec_dev: shell
