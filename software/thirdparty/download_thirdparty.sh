#!/bin/bash


wget http://www.cs.ucl.ac.uk/staff/V.Kolmogorov/software/QPBO-v1.3.src.tar.gz
tar -xzf QPBO-v1.3.src.tar.gz
mv QPBO-v1.3.src/*.cpp QPBO
mv QPBO-v1.3.src/*.h QPBO
mv QPBO-v1.3.src/*.TXT QPBO

wget http://www.f.waseda.jp/hfs/HOCR1.02.zip
unzip HOCR1.02.zip -d tmp
mv tmp/Image.h HOCR/Image.h
mv tmp/HOCR/HOCR.h HOCR/HOCR.h
mv tmp/HOCR/HOCR0.h HOCR/HOCR0.h

rm QPBO-v1.3.src.tar.gz
rm HOCR1.02.zip
rm -r tmp
rm -r QPBO-v1.3.src
