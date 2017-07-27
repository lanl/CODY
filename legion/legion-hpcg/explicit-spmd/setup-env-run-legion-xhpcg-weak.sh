export RXHPCG_EXEC_NAME="legion-xhpcg"

export RXHPCG_START_INDEX=11

export RXHPCG_MAX_SUBBLOCKS=2

export RXHPCG_PPN=1

export RXHPCG_RUN_CMD="aaa -ll:cpu nnn"

export RXHPCG_NUMPE_FUN="X"

export RXHPCG_DATA_DIR_PREFIX="$HOME"

export RXHPCG_NX="16"
export RXHPCG_NY="16"
export RXHPCG_NZ="16"

# Run time in seconds for the benchmark portion of a run.
# Values < 10 result in only one CG set.
export RXHPCG_RT="1"

echo "### run-xhpcg Setup"
env | grep RXHPCG | sort
echo "### run-xhpcg Setup"
