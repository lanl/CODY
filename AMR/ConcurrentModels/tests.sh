#!/bin/bash
   
sh go.sh 60000
sleep 10
echo "DONE WITH GO"

sh d.sh 30000
sleep 10
echo "DONE WITH D"

python postProcess.py -d goTree/ dTree/ -m 32 -o goSpeed1.png goScale1.png goData1.p dSpeed1.png dScale1.png dData1.p


echo "DONE WITH ONE"
sleep 10

sh go.sh 120000
sleep 10
echo "DONE WITH GO"

sh d.sh 60000
sleep 10
echo "DONE WITH D"

python postProcess.py -d goTree/ dTree/ -m 32 -o goSpeed2.png goScale2.png goData2.p dSpeed2.png dScale2.png dData2.p

echo "DONE WITH TWO"
sleep 10

sh go.sh 240000
sleep 10
echo "DONE WITH GO"

sh d.sh 120000
sleep 10
echo "DONE WITH D"

python postProcess.py -d goTree/ dTree/ -m 32 -o goSpeed3.png goScale3.png goData3.p dSpeed3.png dScale3.png dData3.p

echo "DONE WITH THREE"

sleep 10

sh go.sh 480000
sleep 10
echo "DONE WITH GO"

sh d.sh 240000
sleep 10
echo "DONE WITH D"

python postProcess.py -d goTree/ dTree/ -m 32 -o goSpeed4.png goScale4.png goData4.p dSpeed4.png dScale4.png dData4.p

python overhead.py -i goData -o goOverhead
 
 python overhead.py -i dData -o dOverhead















