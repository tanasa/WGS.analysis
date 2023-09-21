#!/bin/bash

FILE=$1

echo "number of deletions :"
grep "DEL" $FILE | wc -l 
echo ""

echo "number of insertions :"
grep "INS" $FILE | wc -l
echo ""

echo "number of inversions :"
grep "INV" $FILE | wc -l
echo ""

echo "number of duplications :"
grep "DUP" $FILE | wc -l
echo ""

echo "number of translocations :"
grep "BND" $FILE | wc -l
echo ""
