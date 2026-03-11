#!/usr/bin/env bash
# =============================================================================
# setup_scratch_volume.sh
#
# Attaches and mounts a dedicated EBS scratch volume for pipeline temp files.
# Run this ONCE per new EC2 instance before starting the pipeline.
#
# Why: /tmp is a 15GB tmpfs (RAM-backed). BAM preprocessing needs ~10-15GB
# per sample. With a scratch EBS volume, multiple samples can be processed
# in parallel without running out of space.
#
# Steps:
#   1. In the AWS Console (or CLI), attach a new gp3 EBS volume to the instance
#      (500GB recommended). It will appear as /dev/xvdf or /dev/nvme1n1.
#   2. Run this script as root or with sudo.
#
# Usage:
#   sudo bash setup_scratch_volume.sh [DEVICE]
#   # DEVICE defaults to auto-detected first unformatted attached disk
# =============================================================================

set -euo pipefail

MOUNT_POINT="/scratch"
DEVICE="${1:-}"

# --- Auto-detect device if not specified -------------------------------------
if [[ -z "${DEVICE}" ]]; then
    # Find first block device that is unformatted (no filesystem)
    for dev in /dev/xvdf /dev/xvdg /dev/nvme1n1 /dev/nvme2n1; do
        if [[ -b "${dev}" ]] && ! blkid "${dev}" &>/dev/null; then
            DEVICE="${dev}"
            break
        fi
    done
fi

if [[ -z "${DEVICE}" ]]; then
    echo "[ERROR] No unformatted EBS volume found. Attach a volume first."
    echo ""
    echo "  AWS CLI: aws ec2 attach-volume --volume-id vol-XXXXXXXXX \\"
    echo "               --instance-id \$(ec2-metadata -i | awk '{print \$2}') \\"
    echo "               --device /dev/xvdf"
    echo ""
    echo "  Or specify device: sudo bash setup_scratch_volume.sh /dev/xvdf"
    exit 1
fi

echo "[setup_scratch] Device   : ${DEVICE}"
echo "[setup_scratch] Mount at : ${MOUNT_POINT}"

# --- Format (only if unformatted) --------------------------------------------
if ! blkid "${DEVICE}" &>/dev/null; then
    echo "[setup_scratch] Formatting ${DEVICE} as ext4..."
    mkfs.ext4 -L scratch "${DEVICE}"
else
    echo "[setup_scratch] ${DEVICE} already formatted — skipping mkfs"
fi

# --- Mount -------------------------------------------------------------------
mkdir -p "${MOUNT_POINT}"

if mountpoint -q "${MOUNT_POINT}"; then
    echo "[setup_scratch] ${MOUNT_POINT} already mounted — skipping"
else
    mount "${DEVICE}" "${MOUNT_POINT}"
    echo "[setup_scratch] Mounted ${DEVICE} at ${MOUNT_POINT}"
fi

# --- Persist across reboots (add to /etc/fstab if not already there) --------
UUID=$(blkid -s UUID -o value "${DEVICE}")
if ! grep -q "${UUID}" /etc/fstab; then
    echo "UUID=${UUID}  ${MOUNT_POINT}  ext4  defaults,nofail  0  2" >> /etc/fstab
    echo "[setup_scratch] Added to /etc/fstab (UUID=${UUID})"
fi

# --- Set permissions ---------------------------------------------------------
chown ec2-user:ec2-user "${MOUNT_POINT}"
chmod 755 "${MOUNT_POINT}"

echo ""
echo "[setup_scratch] Done. Scratch volume ready at ${MOUNT_POINT}."
df -h "${MOUNT_POINT}"
echo ""
echo "  Pipeline TMP_DIR will auto-use /scratch/neoantig_pipeline"
echo "  (config.sh detects /scratch automatically)"
echo ""
echo "  To enable parallel processing, update 01_bam_preprocess.sh:"
echo "    PARALLEL_JOBS=3   # 3 samples concurrently with 500GB scratch"
