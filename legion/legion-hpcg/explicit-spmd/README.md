# A Legion Implementation of HPCG 
Still under development...

## Running
legion-hpcg -ll:cpu [NUMPE] -ll:csize [MEM_IN_B]

## Debugging with Legion Spy
-level legion_spy=2 -logfile log_%.spy

python legion_spy.py -r log_*.spy
python legion_spy.py -lp log_*.spy
python legion_spy.py -ez log_*.spy

-ll:util 0 -hl:separate

-level announce=2
