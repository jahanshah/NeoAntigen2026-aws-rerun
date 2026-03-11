#!/usr/bin/env bash
# =============================================================================
# run_screen.sh — Launch or reattach to a persistent tmux pipeline session
#
# Survives SSH disconnection, AWS session timeouts, and terminal crashes.
#
# Usage:
#   bash run_screen.sh                  # start full pipeline in tmux
#   bash run_screen.sh 2                # start from step 2
#   bash run_screen.sh 3 7              # run steps 3-7
#   bash run_screen.sh attach           # reattach to running session
#   bash run_screen.sh log              # tail the current log (outside tmux)
# =============================================================================

SESSION="neoantig"

case "${1:-}" in

    attach)
        if tmux has-session -t "${SESSION}" 2>/dev/null; then
            echo "[run_screen] Reattaching to tmux session '${SESSION}'..."
            tmux attach-session -t "${SESSION}"
        else
            echo "[run_screen] No session '${SESSION}' found."
            tmux ls 2>/dev/null || echo "No tmux sessions running."
        fi
        exit 0
        ;;

    log)
        LOG="/home/ec2-user/results/bam_preprocessing.log"
        [[ ! -f "${LOG}" ]] && echo "No log found at ${LOG}" && exit 1
        echo "[run_screen] Tailing ${LOG} (Ctrl-C to stop)..."
        tail -f "${LOG}"
        exit 0
        ;;

esac

# --- Check if session already running ----------------------------------------
if tmux has-session -t "${SESSION}" 2>/dev/null; then
    echo "[run_screen] Session '${SESSION}' already running."
    echo "  Reattach : bash run_screen.sh attach"
    echo "  Tail log : bash run_screen.sh log"
    tmux ls
    exit 0
fi

# --- Generate RUN_ID and launch pipeline inside tmux -------------------------
export RUN_ID="res_$(date '+%Y%m%d_%H%M%S')"
LOG="/home/ec2-user/results/bam_preprocessing.log"
mkdir -p /home/ec2-user/results

echo "[run_screen] Starting pipeline in tmux session '${SESSION}'"
echo "  RUN_ID : ${RUN_ID}"
echo "  Log    : ${LOG}"
echo "  Steps  : ${1:-1} → ${2:-12}"
echo ""

tmux new-session -d -s "${SESSION}" -x 220 -y 50 \; \
    send-keys "export RUN_ID=${RUN_ID} && bash /home/ec2-user/code/run_pipeline.sh ${1:-} ${2:-} 2>&1 | tee ${LOG}" Enter

sleep 1
tmux ls
echo ""
echo "  Reattach any time : bash run_screen.sh attach"
echo "  Tail log outside  : bash run_screen.sh log"
echo "  Detach from tmux  : Ctrl-B  D"
