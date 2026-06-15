#!/bin/bash
# run-tests.sh — Run the Allelix test suite in a detached background process.
# Keeps the last 5 runs. No terminal tied up.
#
# Usage:
#   ./scripts/run-tests.sh             # launch tests, returns immediately
#   cat /tmp/allelix-tests/latest.log  # check most recent results
#   ls /tmp/allelix-tests/             # see recent runs

set -euo pipefail

LOG_DIR="/tmp/allelix-tests"
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
KEEP=5

mkdir -p "$LOG_DIR"

# Kill any still-running previous test run
if [ -f "$LOG_DIR/run.pid" ]; then
    OLD_PID=$(cat "$LOG_DIR/run.pid")
    kill "$OLD_PID" 2>/dev/null || true
    rm -f "$LOG_DIR/run.pid"
fi

# Rotate: delete oldest beyond $KEEP (minus 1 to make room for new run)
ls -1t "$LOG_DIR"/run-*.log 2>/dev/null | tail -n +$KEEP | xargs rm -f 2>/dev/null || true

TIMESTAMP=$(date +%Y%m%d-%H%M%S)
LOGFILE="$LOG_DIR/run-$TIMESTAMP.log"

nohup setsid bash -c "
    echo \$\$ > '$LOG_DIR/run.pid'
    cd '$REPO_DIR'
    [ -f '$REPO_DIR/.venv/bin/activate' ] && source '$REPO_DIR/.venv/bin/activate'
    echo '=== Allelix Test Run ===' > '$LOGFILE'
    echo \"Started: \$(date)\" >> '$LOGFILE'
    echo '' >> '$LOGFILE'
    python -m pytest tests/ -x --tb=short --no-cov >> '$LOGFILE' 2>&1
    EXIT_CODE=\$?
    echo '' >> '$LOGFILE'
    echo \"Finished: \$(date)\" >> '$LOGFILE'
    echo \"Exit code: \$EXIT_CODE\" >> '$LOGFILE'
    ln -sf '$LOGFILE' '$LOG_DIR/latest.log'
    rm -f '$LOG_DIR/run.pid'
" > /dev/null 2>&1 &

echo "Tests launched (pid $!). Results → $LOG_DIR/latest.log"
