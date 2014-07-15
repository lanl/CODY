#!/bin/bash
   
sh go.sh 60000
sleep 10
echo "DONE WITH GO"

sh d.sh 30000
sleep 10
echo "DONE WITH D"

python postProcess.py -d goTree/ dStuff/ -m 32 -o goSpeed1.png goOverhead1.png goScale1.png goTime1.pdf goData1.p dSpeed1.png dOverhead1.png dScale1.png dTime1.pdf dData1.p timeData1.p

echo "DONE WITH ONE"
sleep 10

sh go.sh 120000
sleep 10
echo "DONE WITH GO"

sh d.sh 60000
sleep 10
echo "DONE WITH D"

python postProcess.py -d goTree/ dStuff/ -m 32 -o goSpeed2.png goOverhead2.png goScale2.png goTime2.pdf goData2.p dSpeed2.png dOverhead2.png dScale2.png dTime2.pdf dData2.p timeData2.p


echo "DONE WITH TWO"
sleep 10

sh go.sh 240000
sleep 10
echo "DONE WITH GO"

sh d.sh 120000
sleep 10
echo "DONE WITH D"

python postProcess.py -d goTree/ dStuff/ -m 32 -o goSpeed3.png goOverhead3.png goScale3.png goTime3.pdf goData3.p dSpeed3.png dOverhead3.png dScale3.png dTime3.pdf dData3.p timeData3.p


echo "DONE WITH THREE"

sleep 10

sh go.sh 480000
sleep 10
echo "DONE WITH GO"

sh d.sh 240000
sleep 10
echo "DONE WITH D"

python postProcess.py -d goTree/ dStuff/ -m 32 -o goSpeed4.png goOverhead4.png goScale4.png goTime4.pdf goData4.p dSpeed4.png dOverhead4.png dScale4.png dTime4.pdf dData4.p timeData4.p


















