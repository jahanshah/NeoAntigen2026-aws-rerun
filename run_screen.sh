#!/usr/bin/env bash
# =============================================================================
# run_screen.sh — Launch or reattach to a persistent pipeline screen session
#
# Survives SSH disconnection, AWS session timeouts, and terminal crashes.
#
# Usage:
#   bash run_screen.sh                  # start full pipeline in screen
#   bash run_screen.sh 5                # start from step 5
#   bash run_screen.sh 3 7              # run steps 3-7
#   bash run_screen.sh attach           # reattach to running session
#   bash run_screen.sh log              # tail the current log (outside screen)
# =============================================================================

SESSION="neoantig"

case "${1:-}" in

    attach)
        if screen -list | grep -q "${SESSION}"; then
            echo "[run_screen] Reattaching to session '${SESSION}'..."
            screen -r "${SESSION}"
        else
            echo "[run_screen] No session '${SESSION}' found."
            screen -list
        fi
        exit 0
        ;;

    log)
        LOG=$(ls -1t /home/ec2-user/results/bam_preprocessing.log 2>/dev/null | head -1)
        [[ -z "${LOG}" ]] && echo "No log found." && exit 1
        echo "[run_screen] Tailing ${LOG} (Ctrl-C to stop)..."
        tail -f "${LOG}"
        exit 0
        ;;

esac

# --- Check if session already running ----------------------------------------
if screen -list | grep -q "${SESSION}"; then
    echo "[run_screen] Session '${SESSION}' already running."
    echo "  Reattach : bash run_screen.sh attach"
    echo "  Tail log : bash run_screen.sh log"
    screen -list | grep "${SESSION}"
    exit 0
fi

# --- Generate RUN_ID and launch pipeline inside screen -----------------------
export RUN_ID="res_$(date '+%Y%m%d_%H%M%S')"
LOG="/home/ec2-user/results/bam_preprocessing.log"
mkdir -p /home/ec2-user/results

echo "[run_screen] Starting pipeline in screen session '${SESSION}'"
echo "  RUN_ID : ${RUN_ID}"
echo "  Log    : ${LOG}"
echo "  Steps  : ${1:-1} → ${2:-12}"
echo ""
echo "  Reattach any time : bash run_screen.sh attach"
echo "  Tail log outside  : bash run_screen.sh log"
echo ""

screen -dmS "${SESSION}" bash -c "
    export RUN_ID=${RUN_ID}
    bash /home/ec2-user/code/run_pipeline.sh ${1:-} ${2:-} 2>&1 | tee ${LOG}
    echo ''
    echo '[run_screen] Pipeline finished. Press Enter to close this screen.'
    read
"

sleep 1
screen -list | grep "${SESSION}"
echo ""
echo "[run_screen] Detach from screen any time with: Ctrl-A  D"
echo "[run_screen] Reattach with: bash run_screen.sh attach"
